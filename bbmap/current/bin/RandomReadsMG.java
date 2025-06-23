package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import align2.QualityTools;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tax.TaxTree;
import tracker.ReadStats;

/**
 * Generates synthetic reads from one or more input fasta files.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date Feb 8, 2025
 *
 */
public class RandomReadsMG {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		RandomReadsMG x=new RandomReadsMG(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public RandomReadsMG(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 

		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse arguments from the command line */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();

		//Set any necessary Parser defaults here
		//parser.foo=bar;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("in") || a.equals("ref")){
				Tools.getFileOrFiles(b, inputFiles, true, false, false, false);
			}else if(a.equals("depth") || a.equals("cov")){
				minDepth=maxDepth=Float.parseFloat(b);
			}else if(a.equals("mindepth") || a.equals("mincov")){
				minDepth=Float.parseFloat(b);
			}else if(a.equals("maxdepth") || a.equals("maxcov")){
				maxDepth=Float.parseFloat(b);
			}else if(a.equals("depthvariance") || a.equals("variance")){
				depthVariance=Float.parseFloat(b);
			}else if(a.equals("wavecov") || a.equals("sinewave")){
				waveCoverage=Parse.parseBoolean(b);
			}else if(a.equals("waves") || a.equals("sinewaves")){
				numSineWaves=Integer.parseInt(b);
				waveCoverage=numSineWaves>0;
			}else if(a.equals("maxamplitude") || a.equals("waveamp")){
				waveAmp=Float.parseFloat(b);
			}else if(a.equals("oribias")){
				oriBias=Float.parseFloat(b);
			}else if(a.equals("minprob") || a.equals("minwaveprob")){
				minWaveProb=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minPeriod")){
				minPeriod=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("maxPeriod")){
				maxPeriod=Parse.parseIntKMG(b);
				
			}else if(a.equals("mode") || a.equals("depthmode")){
				depthMode=Tools.find(b.toUpperCase(), modes);
				assert(depthMode>=0) : depthMode;
			}else if(a.equals("platform")){
				platform=Tools.find(b.toUpperCase(), platforms);
				assert(platform>=0) : platform;
			}else if(a.equals("insert") || a.equals("avginsert")){
				avgInsert=Float.parseFloat(b);
			}else if(a.equals("paired") || a.equals("int") || a.equals("interleaved")){
				paired=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("reads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("loud")){
				loud=Parse.parseBoolean(b);
			}

			else if(a.equalsIgnoreCase("addErrors")){
				addErrors=Parse.parseBoolean(b);
			}else if(a.equals("qscore") || a.equals("avgq") || a.equals("qavg")){
				meanQScore=Integer.parseInt(b);
			}else if(a.equals("qrange")){
				qScoreRange=Integer.parseInt(b);
			}else if(a.equals("subrate") || a.equals("snprate")){
				subRate=Float.parseFloat(b);
			}else if(a.equals("indelrate")){
				indelRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("pacBioLengthSigma") || a.equals("pbsigma") || a.equals("pacbiosigma")){
				pacBioLengthSigma=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("ontLongTailFactor") || a.equals("tailfactor")){
				ontLongTailFactor=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("srate")){
				sRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("irate")){
				iRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("drate")){
				dRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("hrate")){
				hRate=Float.parseFloat(b);
			}
			
			else if(a.equals("len") || a.equals("length") || a.equals("readlen") || 
					a.equals("readlength") || a.equals("meanlen") || a.equals("avglen")){
				readlen=meanLength=Parse.parseIntKMG(b);
				setReadLength=true;
			}else if(a.equals("maxlen") || a.equals("maxlength")){
				maxLength=Parse.parseIntKMG(b);
				setMaxLength=true;
			}else if(a.equals("minlen") || a.equals("minlength")){
				minLength=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("singleFileThreads") || a.equals("sfthreads") || 
					a.equalsIgnoreCase("sft")){
				singleFileThreads=Integer.parseInt(b);
				Shared.setThreads(Math.max(singleFileThreads, Shared.threads()));
			}else if(a.equalsIgnoreCase("threads") || a.equalsIgnoreCase("t")){
				singleFileThreads=Integer.parseInt(b);
				Shared.setThreads(singleFileThreads);
			}
			
			else if(b!=null && (a.startsWith("cov_") || a.startsWith("depth_"))) {
				float f=Float.parseFloat(b);
				String name=split[0].substring(a.indexOf('_')+1);
				System.err.println("Setting custom depth "+f+" for "+name);
				depthMap.put(name, f);
			}else if(b!=null && Tools.isNumeric(b) && new File(split[0]).isFile()) {
				float f=Float.parseFloat(b);
				inputFiles.add(split[0]);
				String name=ReadWrite.stripPath(split[0]);
				System.err.println("Setting custom depth "+f+" for "+name);
				depthMap.put(name, f);
			}else if(b==null && Tools.find(arg.toUpperCase(), modes)>=0){
				depthMode=Tools.find(arg.toUpperCase(), modes);
				assert(depthMode>=0) : depthMode;
			}else if(b==null && Tools.find(arg.toUpperCase(), platforms)>=0){
				platform=Tools.find(arg.toUpperCase(), platforms);
				assert(platform>=0) : platform;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				File f=new File(arg);
				if(f.exists() && f.canRead()) {
					Tools.getFileOrFiles(arg, inputFiles, true, false, false, false);
				}else {
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
				}
			}
		}

		return parser;
	}

	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}

	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, inputFiles.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		assert(FastaReadInputStream.settingsOK());
	}

	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
		//		assert(false) : "TODO";
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void setPlatformDefaults() {
		if(platform==ILLUMINA) {return;}
		else if(platform==ONT) {
			if(sRate<0) {sRate=0.0025f;}
			if(iRate<0) {iRate=0.0055f;}
			if(dRate<0) {dRate=0.0045f;}
			if(hRate<0) {hRate=0.02f;}
			Shared.FAKE_QUAL=20;
		}else if(platform==PACBIO){
			if(sRate<0) {sRate=0.00015f;}
			if(iRate<0) {iRate=0.000055f;}
			if(dRate<0) {dRate=0.000045f;}
			if(hRate<0) {hRate=0.000015f;}
			Shared.FAKE_QUAL=35;
		}else {
			throw new RuntimeException("Unknown Platform "+platform);
		}
	}
	
	/** Create read streams and process all data */
	void process(Timer t){
		
		setPlatformDefaults();

		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=true;

		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();

		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;

		//Process the reads in separate threads
		spawnThreads(inputFiles, ros);

		if(verbose){outstream.println("Finished; closing streams.");}

		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(ros);

		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;

		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(t.elapsed, readsOut, basesOut, 8));

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Create a Read Input Stream */
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/** Create a Read Output Stream */
	private ConcurrentReadOutputStream makeCros(){
		if(ffout1==null){return null;}

		//Set output buffer size
		final int buff=4;

		//Notify user of output mode
		if(paired && out2==null){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ffout1, ffout2, qfout1, qfout2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}

	/** Spawn process threads */
	private void spawnThreads(final Collection<String> files, final ConcurrentReadOutputStream ros){

		//Do anything necessary prior to processing
		ArrayList<String> flist=new ArrayList<String>(files);
		
		if(singleFileThreads>1 && Shared.threads()>1 && seed<0 && flist.size()==1) {
			//This allows multithreaded processing of a single file.
			//Not really necessary though
			int t=singleFileThreads;
			for(int i=1; i<t; i++) {flist.add(flist.get(0));}
			if(maxReads>0) {maxReads=(maxReads+t-1)/t;}
			else{
				String name=ReadWrite.stripPath(flist.get(0));
				Float depth=depthMap.get(name);
				if(depth==null) {
					depth=randomDepth(Shared.threadLocalRandom(seed))/t;
				}else {
					depth/=t;
				}
				depthMap.put(name, depth);
			}
			loud=false;
		}

		//Determine how many threads may be used
		final int threads=Tools.min(Shared.threads(), flist.size());
//		System.err.println("Using "+threads);
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		AtomicInteger atom=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(flist, ros, i, atom));
		}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}

		//Wait for threads to finish
		waitForThreads(alpt);

		//Do anything necessary after processing

	}

	/** Wait until all worker threads are finished, then return */
	private void waitForThreads(ArrayList<ProcessThread> alpt){

		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			//Accumulate per-thread statistics
			readsProcessed+=pt.readsInT;
			basesProcessed+=pt.basesInT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			success&=pt.success;
		}

		//Track whether any threads failed
		if(!success){errorState=true;}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/**
	 * Determines the coverage depth for a specific input file.
	 * First checks for custom depth settings, then uses the selected distribution mode.
	 * 
	 * @param path Filename to check for custom depth setting
	 * @param taxID Taxonomy ID to check for custom depth setting
	 * @param fnum File number for reporting
	 * @param randy Random number generator
	 * @return The depth to use for this file
	 */
	float chooseDepthForFile(String path, int taxID, int fnum, Random randy) {
		Float custom=null;
		String fname=ReadWrite.stripPath(path);
		if(depthMap.size()>0) {
			custom=depthMap.get(path);
			if(taxID>0 && custom==null) {custom=depthMap.get(Integer.toString(taxID));}
		}
		final float depth;
		if(custom!=null) {depth=custom;}
		else{depth=randomDepth(randy);}
		if(loud) {
			String dstring=(custom==null ? "" : " custom")+
					String.format("depth=%.2f", depth);
			String idstring=taxID>0 ? ("tid "+taxID) : ("name "+fname);
			System.err.println("File "+fnum+", "+idstring+": "+dstring);
		}
		return depth;
	}
	
	float randomDepth(Random randy) {
		final float depth;
		if(depthMode==MIN4) {depth=depthMin4(randy);}
		else if(depthMode==EXP) {depth=depthExp(randy);}
		else if(depthMode==ROOT) {depth=depthRoot(randy);}
		else if(depthMode==LINEAR) {depth=depthLinear(randy);}
		else {throw new RuntimeException("Unknown mode "+depthMode);}
		return depth;
	}

	/**
	 * Generates a depth value using the minimum of 4 random values.
	 * This creates a distribution that is skewed toward lower values.
	 * 
	 * @param randy Random number generator
	 * @return The generated depth value
	 */
	float depthMin4(Random randy) {
		float minRoot=(float)Math.sqrt(minDepth);
		float range=(float)(Math.sqrt(maxDepth)-minRoot);
		final float rootDepth=minRoot+(Tools.min(randy.nextFloat(), randy.nextFloat(), 
				randy.nextFloat(), randy.nextFloat()))*range;
		final float fileDepth=rootDepth*rootDepth;
		return fileDepth;
	}

	/**
	 * Generates a depth value using a uniform linear distribution.
	 * 
	 * @param randy Random number generator
	 * @return The generated depth value
	 */
	float depthLinear(Random randy) {
		float range=(float)(maxDepth-minDepth);
		return randy.nextFloat()*range+minDepth;
	}

	/**
	 * Generates a depth value using a square root distribution.
	 * This creates a distribution that is moderately skewed toward lower values.
	 * 
	 * @param randy Random number generator
	 * @return The generated depth value
	 */
	float depthRoot(Random randy) {
		float range=(float)maxDepth-minDepth;
		float root=randy.nextFloat();
		return root*root*range+minDepth;
	}

	/**
	 * Generates a depth value using an exponential distribution.
	 * This creates a natural long-tailed distribution common in metagenomic samples.
	 * 
	 * @param randy Random number generator
	 * @return The generated depth value
	 */
	float depthExp(Random randy) {
		double lambda=1/Math.sqrt(minDepth*maxDepth);
		double depth=Tools.exponential(randy, lambda);
		while(depth<minDepth || depth>maxDepth) {
			depth=Tools.exponential(randy, lambda);
		}
		return (float)depth;
	}

	public static int mutateIllumina(Read r, int meanQ, int qRange, Random randy) {
		if(r.quality==null) {r.quality=new byte[r.length()];}
		final byte[] bases=r.bases, quals=r.quality;
		int fullRange=qRange*2+1;
		int baseQ=meanQ-qRange;
		int subs=0;
		for(int i=0; i<bases.length; i++) {
			byte b=bases[i];
			if(AminoAcid.isFullyDefined(b)) {
				int q=baseQ+randy.nextInt(fullRange);
				quals[i]=(byte)q;
				float prob=QualityTools.PROB_CORRECT[q];
				if(randy.nextFloat()>prob) {
					int x=AminoAcid.baseToNumber[b];
					x=(x+(randy.nextInt(3)+1))&3;
					bases[i]=AminoAcid.numberToBase[x];
					subs++;
				}
			}else {
				quals[i]=0;
			}
		}
		return subs;
	}
	
	/**
	 * Adds substitutions and deletions to a read at the specified rate, ignoring quality scores.
	 * 
	 * @param r The read to modify
	 * @param rate The probability (0.0-1.0) of introducing a substitution at each position
	 * @param randy Random number generator for randomized decisions
	 * @return The number of indels (insertions + deletions) added to the read
	 */
	public static int addSubs(Read r, float rate, Random randy) {
	    final byte[] bases=r.bases;
	    int subs=0;
	    
	    for(int i=0; i<bases.length; i++) {
	        byte b=bases[i];
	        if(AminoAcid.isFullyDefined(b) && randy.nextFloat()<rate) {
	            // Make a substitution - choose a different base
	            int x=AminoAcid.baseToNumber[b];
	            x=(x+(randy.nextInt(3)+1))&3; // Add 1-3 to avoid same base
	            bases[i]=AminoAcid.numberToBase[x];
	            subs++;
	        }
	    }
	    return subs;
	}

	/**
	 * Adds insertions and deletions to a read at the specified rate, ensuring the final read length
	 * matches the desired length.
	 * <p>
	 * This method introduces indels randomly throughout the read while tracking available "padding"
	 * to ensure the final read is exactly the requested length. Insertions add random bases and
	 * increase available padding. Deletions are only performed when padding is available to prevent
	 * the read from becoming too short. Quality scores for inserted bases are randomly generated
	 * within the specified quality range.
	 * 
	 * @param r The read to modify with indels
	 * @param rate The probability (0.0-1.0) of introducing an indel at each position
	 * @param desiredLength The target length for the modified read
	 * @param meanQ The mean quality score for inserted bases
	 * @param qRange The range around meanQ for quality score randomization
	 * @param randy Random number generator for randomized decisions
	 * @return The number of indels (insertions + deletions) added to the read
	 * 
	 * @throws AssertionError If the resulting read length doesn't match the desired length
	 */
	public static int addIndels(Read r, float rate, int desiredLength, int meanQ, int qRange, Random randy) {
	    final byte[] bases=r.bases;
	    final byte[] quals=r.quality;
	    
	    int padding=r.bases.length-desiredLength;
	    int indels=0;
	    
	    // Create arrays that can accommodate insertions
	    ByteBuilder newBases = new ByteBuilder(desiredLength);
	    ByteBuilder newQuals = quals==null ? null : new ByteBuilder(desiredLength);
	    
		final int fullRange=qRange*2+1;
		final int baseQ=meanQ-qRange;
	    
	    for(int i=0; i<bases.length && newBases.length<desiredLength; i++) {
	        if(randy.nextFloat()<rate) {
	            // Make an indel - 50% insertion, 50% deletion
	            if(randy.nextBoolean()) {
	                // Insertion - add current base plus a random base
	            	newBases.append(bases[i]);
	                if(newQuals!=null) {newQuals.append(quals[i]);}
	                
	                // Insert a random base
	                int x=randy.nextInt(4);
	                newBases.append(AminoAcid.numberToBase[x]);
	                if(newQuals!=null) {
	            		int q=baseQ+randy.nextInt(fullRange);
	                	newQuals.append((byte)q);
	                }
	                indels++;
	                padding++;
	            } else if(padding>0) {
	                // Deletion - skip this base
	            	padding--;
	                indels++;
	                continue; // Skip to next base
	            }
	        } else {
	            // No indel - keep base as is
            	newBases.append(bases[i]);
                if(newQuals!=null) {newQuals.append(quals[i]);}
	        }
	    }
	    
	    // Create right-sized result arrays
	    assert(newBases.length==desiredLength && newBases.array.length==desiredLength);
	    r.bases=newBases.array;
	    if(newQuals!=null) {r.quality=newQuals.array;}
	    return indels;
	}
	
	public int mutateLongRead(Read r, float sRate, float iRate, float dRate, float hRate, Random randy) {
		final float delProb=dRate/Math.max(0.000000000001f, (iRate+dRate));
		final float errProb=sRate+iRate+dRate;
		float bonus=0;
		byte prev=-1;
		int changes=0;
		ByteBuilder bb=new ByteBuilder(r.length()/8+10);
		final byte[] bases=r.bases;
		for(int i=0; i<bases.length; i++) {
			byte b=bases[i];
			if(!AminoAcid.isFullyDefined(b)) {bb.append(b);continue;}

			float f=randy.nextFloat();
			bonus=(b==prev ? bonus+hRate : 0);
			if(f>=errProb+bonus) {
				bb.append(b);
			}else if(f<sRate){//Substitution
				int x=AminoAcid.baseToNumber[b];
				x=(x+(randy.nextInt(3)+1))&3;
				bb.append(AminoAcid.numberToBase[x]);
				changes++;
			}else {//Indel
				if(randy.nextFloat()<delProb) {
					//Deletion, do nothing
				}else {//Insertion
					bb.append(b);
					byte b2=(f>errProb ? b : AminoAcid.numberToBase[randy.nextInt(4)]);
					bb.append(b2);//This is a same-base insertion, sometimes.
				}
				bonus=0;
				changes++;
			}
		}
		r.bases=bb.toBytes();
		return changes;
	}
	
	/**
	 * Generates a random read length following PacBio HiFi distribution.
	 * @param meanLength Mean read length (typically 10000-15000)
	 * @param sigma Standard deviation in log space (typically 0.4-0.5)
	 * @return A randomly generated read length
	 */
	public static int generatePacBioHiFiLength(int minLength, int meanLength, int maxLength, double sigma, Random randy) {
	    // Log-normal distribution
	    double mu = Math.log(meanLength) - 0.5*sigma*sigma; // Adjustment to get desired mean
	    double logLength = mu + randy.nextGaussian()*sigma;
	    int length = (int)Math.round(Math.exp(logLength));
	    
	    // Apply reasonable bounds
	    return Math.max(minLength, Math.min(length, Math.max(meanLength*4, maxLength)));
	}
	
	/**
	 * Generates a random read length following ONT distribution.
	 * @param meanLength Approximate mean read length
	 * @param maxLength Maximum possible read length
	 * @param longTailFactor Controls the heaviness of the tail (0.1-0.3 typical)
	 * @return A randomly generated read length
	 */
	public static int generateONTLength(int minLength, int meanLength, int maxLength, double longTailFactor, Random randy) {
	    // Use mixed distribution approach
	    if(randy.nextDouble() < longTailFactor) {
	        // Generate from heavy tail (exponential)
	        double scale = meanLength * 2; // Scale for longer reads
	        double x = -Math.log(randy.nextDouble()) * scale;
	        return (int)Math.min(x, maxLength);
	    } else {
	        // Generate from log-normal core distribution
	        double sigma = 0.5;
	        double mu = Math.log(meanLength) - 0.5*sigma*sigma;
	        double logLength = mu + randy.nextGaussian()*sigma;
	        return Math.max(minLength, (int)Math.min(Math.exp(logLength), maxLength));
	    }
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private class ProcessThread extends Thread {

		//Constructor
		ProcessThread(final ArrayList<String> files_, final ConcurrentReadOutputStream ros_, 
				int tid_, final AtomicInteger nextFile_){
			files=files_;
			ros=ros_;
			tid=tid_;
			nextFile=nextFile_;
		}

		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			randy=Shared.threadLocalRandom(seed>=0 ? seed+tid : -1);

			//Process the files
			for(int i=nextFile.getAndIncrement(); i<files.size(); i=nextFile.getAndIncrement()) {
				String fname=files.get(i);
				processFile(fname, i);
			}

			//Do anything necessary after processing

			//Indicate successful exit status
			success=true;
		}

		/**
		 * Processes a single input file, generating synthetic reads at the appropriate depth.
		 * 
		 * @param path File name to process
		 * @param fnum File number for tracking and reporting
		 */
		void processFile(String path, int fnum){
//			System.err.println("Thread "+tid+" processing file "+fnum+"; next="+nextFile);
			final String fname=ReadWrite.stripPath(path);
			ConcurrentReadInputStream cris=makeCris(path);

			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			int taxID=TaxTree.parseHeaderStatic2(path, tree);
			if(taxID<0 && ln.size()>0) {
				Read c0=ln.get(0);
				taxID=TaxTree.parseHeaderStatic2(c0.id, tree);
			}
			final float fileDepth=(maxReads>0 ? 0 : chooseDepthForFile(fname, taxID, fnum, randy));
			//			assert(taxID>0) : "Can't parse taxID from "+fname;
			
			if(waveCoverage) {
				covModel=new CoverageModel(numSineWaves, waveAmp, oriBias, minWaveProb, randy, minPeriod, maxPeriod);
			}else {covModel=null;}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
				
				for(Read c : ln) {
					float depthRatio=1f;
					if(varyDepthPerContig) {
						depthRatio=1f+(depthVariance*(randy.nextFloat()+randy.nextFloat()))-depthVariance;
					}
					float contigDepth=depthRatio*fileDepth;
					//					System.err.println("depthRatio = "+depthRatio+"; contigDepth="+contigDepth);
					processContig(c, contigDepth, taxID, fnum, fname);
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				
				//Fetch a new list
				ln=cris.nextList();
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		
		}

		/**
		 * Generates synthetic reads from a single contig at the specified depth.
		 * 
		 * @param contig Source contig to generate reads from
		 * @param depth Coverage depth to generate
		 * @param taxID Taxonomy ID for read headers
		 * @param fnum File number for read headers
		 */
		private void processContig(Read contig, float depth, int taxID, int fnum, String fname) {
			final int basesPerRead=(paired ? 2*readlen : readlen);
			readsInT++;
			basesInT+=contig.length();
			if(contig.length()<basesPerRead+10) {return;}
			if(paired && contig.length()<avgInsert) {return;}

			long basesToGenerate=(maxReads>0 ? maxReads*basesPerRead : (long)(depth*contig.length()));
			long readsGenerated=0;
			long basesGenerated=0;
			ArrayList<Read> list=new ArrayList<Read>(200);
			float variance=varyDepthPerContig ? 0 : randy.nextFloat()*depthVariance;

			//			System.err.println("Generating "+basesToGenerate+" for depth-"+depth+" contig "+contig.id);

			for(long i=0; basesGenerated<basesToGenerate; i++) {
				Read r=generateRead(contig, i, taxID, fnum, contig.numericID, variance, fname);
				if(r!=null) {
					list.add(r);
					readsGenerated+=r.pairCount();
					basesGenerated+=r.pairLength();
				}
				if(list.size()>=200) {
					if(ros!=null) {ros.add(list, 0);}
					list=new ArrayList<Read>(200);
				}
			}
			if(list.size()>0) {if(ros!=null) {ros.add(list, 0);}}
			//			System.err.println("Generated "+basesGenerated+" for depth-"+depth+" contig "+contig.id);

			readsOutT+=readsGenerated;
			basesOutT+=basesGenerated;
		}

		/**
		 * Generates a single or paired read based on the paired flag.
		 * 
		 * @param contig Source contig to generate reads from
		 * @param rnum Read number for identification
		 * @param taxID Taxonomy ID for read headers
		 * @param fnum File number for read headers
		 * @param cnum Contig number for read headers
		 * @param variance Variance value for depth variation
		 * @return The generated read or null if skipped
		 */
		private Read generateRead(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname) {
			if(platform==ILLUMINA) {
				if(paired) {return generateReadPair(contig, rnum, taxID, fnum, cnum, variance, fname);}
				else {return generateReadSingle(contig, rnum, taxID, fnum, cnum, variance, fname);}
			}else if(platform==PACBIO) {return generateLongRead(contig, rnum, taxID, fnum, cnum, variance, fname);}
			else if(platform==ONT) {return generateLongRead(contig, rnum, taxID, fnum, cnum, variance, fname);}
			else {
				throw new RuntimeException("Unknown Platform "+platform);
			}
		}

		/**
		 * Generates a single unpaired read from a contig.
		 * 
		 * @param contig Source contig to generate reads from
		 * @param rnum Read number for identification
		 * @param taxID Taxonomy ID for read headers
		 * @param fnum File number for read headers
		 * @param cnum Contig number for read headers
		 * @param variance Variance value for depth variation
		 * @return The generated read or null if skipped
		 */
		private Read generateLongRead(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname) {
			final int insert;
			if(platform==PACBIO) {
				insert=generatePacBioHiFiLength(minLength, meanLength, maxLength, pacBioLengthSigma, randy);
			}else {
				insert=generateONTLength(minLength, meanLength, maxLength, ontLongTailFactor, randy);
			}
			int paddedLen=insert+(indelRate>0 ? 20 : 0);
			int start=randy.nextInt(contig.length()-paddedLen);
			if(skip((start+paddedLen)/2, contig.length(), variance)
				|| start+paddedLen>contig.length() || start<0){return null;}
			int strand=randy.nextInt(2);
			byte[] bases=Arrays.copyOfRange(contig.bases, start, start+paddedLen);
			if(strand==1) {AminoAcid.reverseComplementBasesInPlace(bases);}
			String header=makeHeader(start, strand, paddedLen, taxID, fnum, cnum, 0, fname);
			Read r=new Read(bases, null, header, rnum);
			if(addErrors) {mutateLongRead(r, sRate, iRate, dRate, hRate, randy);}
			if(subRate>0) {addSubs(r, subRate, randy);}
			if(indelRate>0) {addIndels(r, indelRate, paddedLen, meanQScore, qScoreRange, randy);}
			return r;
		}

		/**
		 * Generates a single unpaired read from a contig.
		 * 
		 * @param contig Source contig to generate reads from
		 * @param rnum Read number for identification
		 * @param taxID Taxonomy ID for read headers
		 * @param fnum File number for read headers
		 * @param cnum Contig number for read headers
		 * @param variance Variance value for depth variation
		 * @return The generated read or null if skipped
		 */
		private Read generateReadSingle(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname) {
			final int paddedLen=readlen+(indelRate>0 ? 5 : 0);
			int insert=paddedLen;
			int start=randy.nextInt(contig.length()-insert);
			if(skip((start+insert)/2, contig.length(), variance) 
					|| start+paddedLen>contig.length() || start<0){return null;}
			int strand=randy.nextInt(2);
			byte[] bases=Arrays.copyOfRange(contig.bases, start, start+paddedLen);
			if(strand==1) {AminoAcid.reverseComplementBasesInPlace(bases);}
			String header=makeHeader(start, strand, insert, taxID, fnum, cnum, 0, fname);
			Read r=new Read(bases, null, header, rnum);
			if(addErrors) {mutateIllumina(r, meanQScore, qScoreRange, randy);}
			if(subRate>0) {addSubs(r, subRate, randy);}
			if(indelRate>0) {addIndels(r, indelRate, readlen, meanQScore, qScoreRange, randy);}
			return r;
		}

		/**
		 * Generates a paired read from a contig with appropriate insert size.
		 * 
		 * @param contig Source contig to generate reads from
		 * @param rnum Read number for identification
		 * @param taxID Taxonomy ID for read headers
		 * @param fnum File number for read headers
		 * @param cnum Contig number for read headers
		 * @param variance Variance value for depth variation
		 * @return The generated read pair or null if skipped
		 */
		private Read generateReadPair(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname) {
			double g=randy.nextGaussian()*0.25f;
			int insert=(int)((1+g)*avgInsert);
			final int paddedLen=readlen+(indelRate>0 ? 5 : 0);
			while(insert>=contig.length() || insert<paddedLen) {
				g=randy.nextGaussian()*0.25f;
				insert=(int)((1+g)*avgInsert);
			}
			int start1=randy.nextInt(contig.length()-insert);
			int strand=randy.nextInt(2);
			int start2=start1+insert-paddedLen;
			if(skip((start1+insert)/2, contig.length(), variance)
				|| start1+paddedLen>contig.length() || start1<0
				|| start2+paddedLen>contig.length() || start2<0){return null;}

			byte[] bases1=Arrays.copyOfRange(contig.bases, start1, start1+paddedLen);
			byte[] bases2=Arrays.copyOfRange(contig.bases, start2, start2+paddedLen);
			AminoAcid.reverseComplementBasesInPlace(bases2);
			if(strand==1) {
				byte[] temp=bases1;
				bases1=bases2;
				bases2=temp;
			}
			String header1=makeHeader(start1, strand, insert, taxID, fnum, cnum, 0, fname);
			String header2=makeHeader(start1, strand, insert, taxID, fnum, cnum, 1, fname);
			Read r1=new Read(bases1, null, header1, rnum);
			Read r2=new Read(bases2, null, header2, rnum);
			r2.setPairnum(1);
			r1.mate=r2;
			r2.mate=r1;
			if(addErrors) {
				mutateIllumina(r1, meanQScore, qScoreRange, randy);
				mutateIllumina(r2, meanQScore, qScoreRange, randy);
			}
			if(subRate>0) {
				addSubs(r1, subRate, randy);
				addSubs(r2, subRate, randy);
			}
			if(indelRate>0) {
				addIndels(r1, indelRate, readlen, meanQScore, qScoreRange, randy);
				addIndels(r2, indelRate, readlen, meanQScore, qScoreRange, randy);
			}
			return r1;
		}

		/**
		 * Determines whether to skip generating a read at a particular position.
		 * This is used to create within-contig coverage variation.
		 * 
		 * @param midpoint Position in the contig
		 * @param clen Length of the contig
		 * @param variance Variance parameter controlling the skip probability
		 * @return True if this position should be skipped
		 */
		private boolean skip(int midpoint, int clen, float variance) {
			if(covModel!=null) {
				return !covModel.shouldGenerateReadAt(midpoint, clen, randy);
			}else if(variance<=0) {return false;}
			
			//Linear model only
			float maxSkipProb=1f-1f/(1f+variance);
			float relativePosition=midpoint/(float)clen;
			float skipProb=relativePosition*maxSkipProb;
			return randy.nextFloat()<skipProb;
		}

		/**
		 * Creates a read header with genome source information encoded.
		 * Format: f_[filenum]_c_[contignum]_s_[strand]_p_[position]_i_[insert]_tid_[taxid] [pairnum]:
		 * 
		 * @param start Starting position in the contig
		 * @param strand Read strand (0=forward, 1=reverse)
		 * @param insert Insert size for paired reads
		 * @param taxID Taxonomy ID
		 * @param fnum File number
		 * @param cnum Contig number
		 * @param pnum Pair number (0=first, 1=second)
		 * @return Formatted header string
		 */
		private String makeHeader(int start, int strand, int insert, int taxID, 
				int fnum, long cnum, int pnum, String fname) {
			bb.clear().append('f').under().append(fnum).under().append('c').under().append(cnum);
			bb.under().append('s').under().append(strand).under().append('p').under().append(start);
			bb.under().append('i').under().append(insert).under();
			if(taxID>0) {bb.append("tid").under().append(taxID);}
			else {bb.append("name").under().append(fname);}
			bb.space().append(pnum+1).colon();
			return bb.toString();
		}

		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;

		/** Number of input contigs processed by this thread */
		protected long readsInT=0;
		/** Number of input bases processed by this thread */
		protected long basesInT=0;

		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
		final AtomicInteger nextFile;
		private final ArrayList<String> files;
		private Random randy;
		private CoverageModel covModel;

		private ByteBuilder bb=new ByteBuilder(128);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private LinkedHashSet<String> inputFiles=new LinkedHashSet<String>();

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;

	/** Override output file extension */
	private String extout=null;

	private String taxTreeFile=null;

	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;
	
	private boolean loud=true;

	private AtomicLong nextReadID=new AtomicLong(0);

	private float minDepth=1;
	private float maxDepth=256;
	private float depthVariance=0.5f;
	
	private boolean waveCoverage=false;
	private int numSineWaves=4;
	private float waveAmp=0.7f;
	private float oriBias=0.25f;
	private float minWaveProb=0.1f;
	private int minPeriod=2000;
	private int maxPeriod=80000;
	
	private long seed=-1;
	private boolean varyDepthPerContig=false;
	private HashMap<String, Float> depthMap=new HashMap<String, Float>();
	private long maxReads=-1;
	private float subRate=0;
	private float indelRate=0;

	private float avgInsert=300;
	private int readlen=150;
	private boolean paired=true;
	private boolean addErrors=false;
	private int meanQScore=25;
	private int qScoreRange=0;
	
	private float pacBioLengthSigma=0.5f;
	private float ontLongTailFactor=0.2f;
	private int minLength=1000;
	private int meanLength=15000;
	private int maxLength=100000;
	private int singleFileThreads=1;
	
	private float sRate=-1;
	private float iRate=-1;
	private float dRate=-1;
	private float hRate=-1;

	static final String[] modes={"MIN4", "EXP", "ROOT", "LINEAR"};
	static final int MIN4=0, EXP=1, ROOT=2, LINEAR=3;
	int depthMode=MIN4;

	static final String[] platforms={"ILLUMINA", "ONT", "PACBIO"};
	static final int ILLUMINA=0, ONT=1, PACBIO=2;
	int platform=ILLUMINA;
	
	boolean setReadLength=false;
	boolean setMaxLength=false;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;

	private final TaxTree tree;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;

}
