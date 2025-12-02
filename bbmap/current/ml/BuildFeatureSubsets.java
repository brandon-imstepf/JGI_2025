package ml;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.LineParser1;
import structures.ByteBuilder;

/**
 * ByteStream/LineParser implementation of the feature-subset slicer.
 * Streams the gzipped .all/.positives/.negatives TSVs for a dataset and writes the
 * subset variants (full, floats-only, onehot-only, floats+outer-onehot) without ever
 * touching String.split().
 */
public class BuildFeatureSubsets {

	public static void main(String[] args){
		Timer t=new Timer();
		BuildFeatureSubsets bfs=new BuildFeatureSubsets(args);
		bfs.process(t);
		Shared.closeStream(bfs.outstream);
	}

	public BuildFeatureSubsets(String[] args){
		PreParser pp=new PreParser(args, null, false);
		args=pp.args;
		outstream=pp.outstream;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		final Parser parser=parse(args);
		overwrite=parser.overwrite;
		append=parser.append;
		validateParams();
	}

	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0;i<args.length;i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("in") || a.equals("indir")){
				inDir=b;
			}else if(a.equals("out") || a.equals("outdir")){
				outDir=b;
			}else if(a.equals("subsets")){
				if(b!=null){
					subsetNames=b.split(",");
				}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//handled by parser
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		return parser;
	}

	private void validateParams(){
		if(inDir==null){throw new RuntimeException("Input directory required: in=<dir>");}
		if(outDir==null){throw new RuntimeException("Output directory required: out=<dir>");}
		inputDirectory=new File(inDir);
		if(!inputDirectory.isDirectory()){
			throw new RuntimeException("Input directory "+inputDirectory+" does not exist");
		}
		outputDirectory=new File(outDir);
		if(!outputDirectory.exists() && !outputDirectory.mkdirs()){
			throw new RuntimeException("Unable to create output directory "+outputDirectory);
		}
	}

	private void process(Timer t){
		File allFile=findSingle(".all.tsv");
		File posFile=findSingle(".positives.tsv");
		File negFile=findSingle(".negatives.tsv");
		List<File> sources=Arrays.asList(allFile, posFile, negFile);
		if(verbose){
			outstream.println("Input files:");
			for(File f : sources){outstream.println("  "+f.getAbsolutePath());}
		}
		int totalCols=inferTotalColumns(allFile);
		labelIndex=totalCols-1;
		int[][] masks=computeMasks(totalCols, subsetNames);
		for(File src : sources){
			processFile(src, masks);
		}
		t.stop();
		outstream.println("Finished in "+t);
	}

	private void processFile(File src, int[][] masks){
		if(src==null){return;}
		ByteFile bf=ByteFile.makeByteFile(src.getAbsolutePath(), true);
		LineParser1 parser=new LineParser1('\t');
		ByteStreamWriter[] writers=openWriters(src.getName(), masks.length);
		ByteBuilder[] builders=new ByteBuilder[masks.length];
		for(int i=0;i<builders.length;i++){builders[i]=new ByteBuilder(512);}
		byte[] line=bf.nextLine();
		while(line!=null){
			final int effectiveLen=trimLine(line);
			if(effectiveLen<=0){
				line=bf.nextLine();
				continue;
			}
			parser.setLimited(line, effectiveLen);
			if(parser.terms()<=labelIndex){
				line=bf.nextLine();
				continue;
			}
			final byte label=sanitizeLabel(line, parser, labelIndex);
			if(label!='0' && label!='1'){
				line=bf.nextLine();
				continue;
			}
			for(int s=0;s<masks.length;s++){
				ByteBuilder bb=builders[s];
				bb.clear();
				int[] mask=masks[s];
				for(int j=0;j<mask.length;j++){
					if(j>0){bb.append('\t');}
					int col=mask[j];
					if(col==labelIndex){
						bb.append(label);
					}else{
						appendTerm(bb, line, parser, col);
					}
				}
				bb.append('\n');
				writers[s].print(bb);
			}
			line=bf.nextLine();
		}
		errorState|=bf.close();
		for(ByteStreamWriter bsw : writers){
			if(bsw!=null){errorState|=bsw.poisonAndWait();}
		}
	}

	private ByteStreamWriter[] openWriters(String baseName, int count){
		ByteStreamWriter[] arr=new ByteStreamWriter[count];
		String suffix=baseName.endsWith(".tsv.gz") ? ".tsv.gz" : baseName.endsWith(".tsv") ? ".tsv" : "";
		for(int i=0;i<count;i++){
			String subset=subsetNames[i];
			File dir=new File(outputDirectory, subset);
			if(!dir.exists()){dir.mkdirs();}
			String outName;
			if(suffix.length()>0 && baseName.endsWith(suffix)){
				outName=baseName.substring(0, baseName.length()-suffix.length())+".subset_"+subset+suffix;
			}else{
				outName=baseName+".subset_"+subset+".tsv.gz";
			}
			File outFile=new File(dir, outName);
			FileFormat ff=FileFormat.testOutput(outFile.getAbsolutePath(), FileFormat.TXT, null, true, overwrite, append, false);
			arr[i]=ByteStreamWriter.makeBSW(ff);
		}
		return arr;
	}

	private int inferTotalColumns(File file){
		ByteFile bf=ByteFile.makeByteFile(file.getAbsolutePath(), true);
		byte[] line=bf.nextLine();
		while(line!=null && trimLine(line)<=0){
			line=bf.nextLine();
		}
		if(line==null){
			bf.close();
			throw new RuntimeException("Could not read from "+file.getAbsolutePath());
		}
		final int len=trimLine(line);
		int cols=1;
		for(int i=0;i<len;i++){
			if(line[i]=='\t'){cols++;}
		}
		bf.close();
		return cols;
	}

	private static int[][] computeMasks(int totalCols, String[] subsets){
		int labelIdx=totalCols-1;
		int onehotTotal=Math.max(0, labelIdx-FLOAT_COUNT);
		int segment=onehotTotal/3;
		int remainder=onehotTotal-segment*3;
		int startBlock=FLOAT_COUNT;
		int firstEnd=startBlock+segment;
		int lastStart=startBlock+segment+remainder+segment;
		int[][] masks=new int[subsets.length][];
		for(int i=0;i<subsets.length;i++){
			String name=subsets[i];
			ArrayList<Integer> idxs=new ArrayList<>();
			if("full".equals(name)){
				for(int c=0;c<totalCols;c++){idxs.add(c);}
			}else if("floats_only".equals(name)){
				for(int c=0;c<FLOAT_COUNT;c++){idxs.add(c);}
				idxs.add(labelIdx);
			}else if("onehot_only".equals(name)){
				for(int c=startBlock;c<labelIdx;c++){idxs.add(c);}
				idxs.add(labelIdx);
			}else if("floats_plus_onehot_no_middle".equals(name)){
				for(int c=0;c<FLOAT_COUNT;c++){idxs.add(c);}
				for(int c=startBlock;c<firstEnd;c++){idxs.add(c);}
				for(int c=lastStart;c<labelIdx;c++){idxs.add(c);}
				idxs.add(labelIdx);
			}else{
				for(int c=0;c<totalCols;c++){idxs.add(c);}
			}
			int[] mask=new int[idxs.size()];
			for(int j=0;j<mask.length;j++){mask[j]=idxs.get(j);}
			masks[i]=mask;
		}
		return masks;
	}

	private static void appendTerm(ByteBuilder bb, byte[] line, LineParser1 parser, int term){
		parser.setBounds(term);
		final int start=parser.a();
		final int end=parser.b();
		for(int i=start;i<end;i++){
			byte c=line[i];
			if(c==0 || c=='\r' || c=='\n'){continue;}
			bb.append(c);
		}
	}

	private static byte sanitizeLabel(byte[] line, LineParser1 parser, int term){
		parser.setBounds(term);
		final int start=parser.a();
		final int end=parser.b();
		byte result='?';
		for(int i=start;i<end;i++){
			byte c=line[i];
			if(c=='0' || c=='1'){result=c;}
		}
		return result;
	}

	private static int trimLine(byte[] line){
		int len=line.length;
		while(len>0){
			byte b=line[len-1];
			if(b==0 || b=='\n' || b=='\r'){
				len--;
			}else{
				break;
			}
		}
		return len;
	}

	private File findSingle(String token){
		File[] matches=inputDirectory.listFiles((dir, name) -> name.endsWith(".tsv.gz") && name.contains(token));
		if(matches==null || matches.length==0){
			throw new RuntimeException("Could not find file containing '"+token+"' in "+inputDirectory);
		}
		Arrays.sort(matches);
		if(matches.length>1 && verbose){
			outstream.println("Warning: multiple files matched "+token+"; using "+matches[0].getName());
		}
		return matches[0];
	}

	private static final int FLOAT_COUNT=8;
	private static final String[] DEFAULT_SUBSETS=new String[]{"full","floats_only","onehot_only","floats_plus_onehot_no_middle"};

	private String[] subsetNames=DEFAULT_SUBSETS.clone();
	private String inDir=null;
	private String outDir=null;
	private File inputDirectory=null;
	private File outputDirectory=null;
	private int labelIndex=-1;

	private PrintStream outstream=System.err;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean verbose=false;
	public boolean errorState=false;
}
