package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.locks.ReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.GffLine;
import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import prok.ProkObject;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.FloatList;
import structures.IntHashMap;
import structures.IntLongHashMap;
import structures.ListNum;
import structures.LongList;
import tax.Lineage;
import tax.TaxNode;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Evaluates the quality of metagenomic bins based on various metrics including
 * completeness, contamination, taxonomic assignment, and gene content.
 * <p>
 * This tool can process multiple bin files (FASTA format) and generate comprehensive
 * quality reports including completeness/contamination metrics, taxonomic classification,
 * and gene annotation. It supports integration with external tools like CheckM, EukCC,
 * GTDB-Tk, and CAMI for validation, and can optionally identify genes and RNA sequences.
 * <p>
 * Quality categories include:
 * - UHQ (Ultra-High Quality): >99% completeness, <1% contamination
 * - VHQ (Very High Quality): >95% completeness, <2% contamination
 * - HQ (High Quality): >90% completeness, <5% contamination (with RNA genes if useRNA=true)
 * - MQ (Medium Quality): >50% completeness, <10% contamination
 * - LQ (Low Quality): Lower quality bins not meeting MQ criteria
 * - VLQ (Very Low Quality): <20% completeness or >20% contamination
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date Feb 8, 2025
 */
public class GradeBins implements Accumulator<GradeBins.ProcessThread> {
	
	/**
	 * Entry point for the GradeBins program.
	 * Initializes the program, processes bins, and reports results.
	 * 
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		GradeBins x=new GradeBins(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Creates a new GradeBins instance and parses command-line arguments.
	 * Loads necessary reference data and initializes processing parameters.
	 * 
	 * @param args Command-line arguments
	 */
	public GradeBins(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("in")){
				Tools.getFileOrFiles(b, in, true, false, false, false);
			}else if(a.equals("size")){
				totalSize=Parse.parseKMG(b);
			}else if(a.equals("minsize")){
				minSize=Parse.parseIntKMG(b);
			}else if(a.equals("ref") || a.equals("contigs") || a.equals("assembly")){
				ref=b;
			}else if(a.equals("hist")){
				hist=b;
			}else if(a.equalsIgnoreCase("contamHist")){
				contamHist=b;
			}else if(a.equals("ccplot")){
				ccplot=b;
			}
			
			else if(a.equals("report")){
				report=b;
			}else if(a.equals("taxin")){
				taxIn=b;
			}else if(a.equals("taxout")){
				taxOut=b;
			}else if(a.equals("tax") || a.equals("size")){
				tax=b;
			}else if(a.equals("cov")){
				cov=b;
			}else if(a.equals("loadmt")){
				loadMT=Parse.parseBoolean(b);
			}else if(a.equals("tree") || a.equals("usetree")){
				if(b==null || Parse.isBoolean(b)) {useTree=Parse.parseBoolean(b);}
				else if(new File(b).exists()) {
					BinObject.treePath=b;
					useTree=true;
				}else {
					assert(false) : "Bad argument: "+arg;
				}
			}
			
			else if(a.equalsIgnoreCase("checkm")){
				checkMFile=b;
			}else if(a.equalsIgnoreCase("eukcc")){
				eukCCFile=b;
			}else if(a.equalsIgnoreCase("cami")){
				camiFile=b;
			}else if(a.equalsIgnoreCase("gtdb") || a.equalsIgnoreCase("gtdbtk")){
				gtdbFile=b;
			}else if(a.equalsIgnoreCase("pgm")){
				GeneTools.pgmFile=b;
			}else if(a.equalsIgnoreCase("gff")){
				gffFile=b;
			}else if(a.equalsIgnoreCase("imgmap")){
				imgMapFile=b;
			}else if(a.equalsIgnoreCase("spectra")){
				spectraFile=b;
			}else if(a.equalsIgnoreCase("quickclade")){
				runQuickClade=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("callgenes")){
				callGenes=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("userna") || a.equals("rna") || a.equals("ribo")){
				useRNA=Parse.parseBoolean(b);
			}else if(a.equals("aligner") || a.equals("idaligner")){
				GeneCaller.useIDAligner=(b==null || !("f".equals(b) || "false".equals(b)));
				if(GeneCaller.useIDAligner) {aligner.Factory.setType(b);}
			}else if(b==null && new File(arg).isFile()){
//				FileFormat.PRINT_WARNING=false;
				FileFormat ff=FileFormat.testInput(arg, FileFormat.TXT, null, false, false);
//				FileFormat.PRINT_WARNING=true;
				String lc=arg.toLowerCase();
				if(ff.fasta()) {
					in.add(arg);
				}else if(ff.pgm()) {
					GeneTools.pgmFile=arg;
				}else if(ff.gff()) {
					gffFile=arg;
				}else if(ff.clade()) {
					spectraFile=arg;
				}else if(lc.contains("checkm") && checkMFile==null) {
					checkMFile=arg;
				}else if(lc.contains("cami") && camiFile==null) {
					camiFile=arg;
				}else if(lc.contains("gtdb") && gtdbFile==null) {
					gtdbFile=arg;
				}else if(lc.contains("eukcc") && eukCCFile==null) {
					eukCCFile=arg;
				}else if(lc.equals("tax.txt") && taxIn==null) {
//					System.err.println("Adding tax "+arg);
					taxIn=arg;
				}else if(DataLoader.looksLikeCovFile(arg) && cov==null) {
//					System.err.println("Adding cov "+arg);
					cov=arg;
				}else {
//					System.err.println("Adding bin "+arg);
					in.add(arg);
				}
			}else if(b==null && new File(arg).isDirectory()){
				Tools.getFileOrFiles(arg, in, true, false, false, false);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
//			out1=parser.out1;
		}
		
		if(callGenes) {
			GeneTools.loadPGM();
			CallGenes.callCDS=CallGenes.calltRNA=CallGenes.call16S=
					CallGenes.call23S=CallGenes.call5S=CallGenes.call18S=true;
		}
		loadGff();
		loadSpectra();
		loadCov();
		makeLevelMaps();
		if(gtdbFile!=null || (cladeIndex!=null && report!=null)) {useTree=true;}
		if(useTree) {BinObject.loadTree();}
	}
	
	/**
	 * Initializes taxonomy level maps for tracking classified organisms.
	 * Creates maps for all bins, medium-quality bins, and high-quality bins.
	 */
	static synchronized void makeLevelMaps() {
		if(levelMaps!=null) {return;}
		levelMaps=new IntHashMap[TaxTree.LIFE+1];
		levelMapsMQ=new IntHashMap[TaxTree.LIFE+1];
		levelMapsHQ=new IntHashMap[TaxTree.LIFE+1];
		for(int i=0; i<levelMaps.length; i++) {
			levelMaps[i]=new IntHashMap();
			levelMapsMQ[i]=new IntHashMap();
			levelMapsHQ[i]=new IntHashMap();
		}
	}
	
	/**
	 * Loads coverage information for contigs.
	 * Called once to initialize the coverage map.
	 */
	static synchronized void loadCov() {
		if(cov==null || covMap!=null) {return;}
		covMap=DataLoader.loadCovFile(cov);
	}
	
	/**
	 * Loads gene annotations from GFF files.
	 * Called once to initialize the GFF map.
	 */
	static synchronized void loadGff() {
		if(gffFile==null || gffMap!=null) {return;}
		HashMap<String, String> imgMap=loadImgMap(imgMapFile);
		System.err.println("Loading "+gffFile);
		ArrayList<GffLine> lines=GffLine.loadGffFile(gffFile, "rRNA,tRNA", callGenes);
		gffMap=new HashMap<String, ArrayList<GffLine>>();
		for(GffLine line : lines) {
			if(imgMap!=null) {
				String key=line.seqid;
				String value=imgMap.get(key);
				if(value!=null) {line.seqid=value;}
			}
			ArrayList<GffLine> value=gffMap.get(line.seqid());
			if(value==null) {gffMap.put(line.seqid(), value=new ArrayList<GffLine>(2));}
			value.add(line);
		}
	}
	
	/**
	 * Loads mapping between IMG IDs and contig names.
	 * 
	 * @param fname IMG mapping file
	 * @return Mapping between IMG IDs and contig names
	 */
	static HashMap<String, String> loadImgMap(String fname){
		if(fname==null) {return null;}
		HashMap<String,String> map=new HashMap<String,String>();
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		LineParser1 lp=new LineParser1('\t');
		for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()) {
			for(byte[] line : ln) {
				lp.set(line);
				String a=lp.parseString(0);
				String b=lp.parseString(1);
				String old=map.put(a, b);
				String old2=map.put(b, a);
				assert(old==null) : "Evicted "+old+" for "+a+" -> "+b;
			}
		}
		return map;
	}
	
	/**
	 * Loads the k-mer composition database for taxonomy assignment.
	 * Initializes the QuickClade index if needed.
	 */
	static void loadSpectra() {
		if(runQuickClade && spectraFile==null) {spectraFile=CladeSearcher.defaultRef();}
		if(spectraFile!=null) {runQuickClade=true;}
		if(spectraFile==null || cladeIndex!=null) {return;}
		cladeIndex=CladeIndex.loadIndex(spectraFile);
	}
	
	/**
	 * Main processing method that loads bins, calculates quality metrics,
	 * and generates reports.
	 * 
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		BinObject.grading=true;
		
		if(tax!=null && taxIn==null && taxOut==null) {
			boolean taxExists=(tax==null ? false : new File(tax).canRead());
			if(taxExists && ref==null) {taxIn=tax;}
			else {taxOut=tax;}
//			assert(false) : taxExists+", "+ref+", "+taxIn;
		}
		
		if(taxIn!=null) {
			if(ref!=null) {
				System.err.println("Reading from "+taxIn+" instead of "+ref);
			}
			sizeMap=loadTaxIn(taxIn);
		}else {
			sizeMap=makeSizeMap(ref);
		}
		System.err.println("Made size map.");
		if(taxOut!=null) {
			writeTaxOut(taxOut, sizeMap, countMap);
		}
		checkMMap=loadCheckM(checkMFile);
		eukCCMap=loadEukCC(eukCCFile);
		camiMap=loadCami(camiFile);
		gtdbMap=loadGTDBDir(gtdbFile);
		Timer t2=new Timer(System.err, false);
		System.err.print("Loading bins: ");
		ArrayList<BinStats> bins=(loadMT ? loadMT(in) : loadST(in));
		t2.stopAndPrint();
		
		printResults(bins);
		
		t.stop();
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
	}
	
	/**
	 * Adds taxonomic information to a bin at different taxonomic levels.
	 * Updates taxonomy counters for reporting.
	 * 
	 * @param bin BinStats object to annotate
	 */
	private void addTaxLevels(BinStats bin) {
		Lineage lineage=null;
		if(gtdbMap!=null) {
			lineage=gtdbMap.get(bin.name);
			if(lineage!=null) {
				bin.lineage=lineage.line;
			}
		}
		int tid=bin.taxid;
		if(tid<1) {tid=bin.taxid=TaxTree.LIFE_ID;}
		if(lineage==null) {lineage=new Lineage(tid);}
		addTaxLevels(bin, lineage);
	}
	
	/**
	 * Adds taxonomic information to a bin using a specific lineage.
	 * Updates taxonomy counters for reporting.
	 * 
	 * @param bin BinStats object to annotate
	 * @param lineage Taxonomic lineage to assign
	 */
	private void addTaxLevels(BinStats bin, Lineage lineage) {
		boolean hq=bin.hq(useRNA);
		boolean mq=bin.mq(useRNA);
//		System.err.println("Incrementing lineage for "+bin.taxid);
		for(TaxNode node : lineage.nodes) {
			if(node!=null) {
//				System.err.print('.');
				levelMaps[node.level].increment(node.id);
				if(mq) {levelMapsMQ[node.level].increment(node.id);}
				if(hq) {levelMapsHQ[node.level].increment(node.id);}
			}
		}
//		assert(false);
	}
	
	/**
	 * Prints overall results for the processed bins.
	 * Includes quality metrics, taxonomic distribution, and score calculations.
	 * 
	 * @param bins List of BinStats to report on
	 */
	void printResults(ArrayList<BinStats> bins) {
		for(BinStats bin : bins) {
			readsProcessed+=bin.contigs;
			basesProcessed+=bin.size;
			sizes.add(bin.size);
			if(useTree) {addTaxLevels(bin);}
		}
		
		if(verbose){outstream.println("Finished.");}

		outstream.println();
		printCleanDirty(bins);
		
		outstream.println();
		printL90(sizes, totalSize);

		outstream.println();
		printScore(bins, totalSize, totalContigs, taxIDsIn, true);
		
		outstream.println();
		printBinQuality(bins, minSize, useRNA, outstream);
		
		if(useTree) {
			outstream.println();
			printTaxLevels(bins, outstream);
		}
		
		if(hist!=null) {
			ChartMaker.makeChartFromBinStats(hist, bins);
		}
		if(ccplot!=null) {
			ChartMaker.writeCCPlot(ccplot, bins);
		}
		if(contamHist!=null) {
			ChartMaker.writeContamHist(contamHist, bins);
		}
		if(report!=null) {
			printClusterReport(bins, minSize, report);
		}
	}
	
	/**
	 * Calculates and prints bin score metrics including completeness,
	 * contamination, and overall quality score.
	 * 
	 * @param bins List of BinStats to score
	 * @param totalSize Total size of the reference dataset
	 * @param totalContigs Total number of contigs in the reference
	 * @param taxIDsIn Number of unique taxonomy IDs in the reference
	 * @param validation Whether to print validation metrics
	 */
	public static void printScore(ArrayList<BinStats> bins, 
			long totalSize, long totalContigs, long taxIDsIn, boolean validation) {
		long cleanContigs=0, contamContigs=0;
		long cleanSize=0, contamSize=0;
		long badContigs=0;
		double compltScore=0, contamScore=0;
		double totalScore=0, totalScore2=0;
		IntHashMap tidBins=new IntHashMap();
		int labels=0;
		for(BinStats bin : bins) {
			if(bin.taxid>0) {
				tidBins.increment(bin.taxid);
				labels++;
			}
			long contam=Math.round(bin.contam*bin.size);
			contamScore+=contam;
			compltScore+=Math.round(bin.complt*(bin.size-contam));
			double score=Math.max(0, bin.complt-5*bin.contam);
			totalScore+=score;
			totalScore2+=score*score;
			badContigs+=bin.badContigs;
			if(contam<1) {
				cleanSize+=bin.size;
				cleanContigs+=bin.contigs;
			}else {
				contamSize+=bin.size;
				contamContigs+=bin.contigs;
			}
		}

		outstream.println("Sequence Recovery:           \t"+
				String.format("%.3f", (cleanSize+contamSize)*100.0/totalSize));
		outstream.println("Contig Recovery:             \t"+
				String.format("%.3f", (cleanContigs+contamContigs)*100.0/totalContigs));
		if(taxIDsIn>0) {
			outstream.println("Bad Contigs:                 \t"+
					String.format("%.3f", badContigs*100.0/(cleanContigs+contamContigs)));
			outstream.println("Genomes Represented:         \t"+
					String.format("%.3f", (tidBins.size())*100.0/taxIDsIn));
		}
		if(validation) {
			outstream.println("Completeness Score:          \t"+
					String.format("%.3f", 100*compltScore/totalSize));
			outstream.println("Contamination Score:         \t"+
					String.format("%.4f", 100*contamScore/totalSize));
//			outstream.println("Total Score:                 \t"+
//					String.format("%.2f", totalScore));
			outstream.println("Total Score:               \t"+
					String.format("%.2f", totalScore2));
		}
	}
	
	/**
	 * Generates a score string for a collection of bins.
	 * 
	 * @param bins Collection of bins to score
	 * @param minSize Minimum size threshold for bins to include
	 * @param sizeMap Size map for contamination calculation
	 * @return Formatted score string
	 */
	static String toScoreString(ArrayList<? extends Bin> bins, int minSize, IntLongHashMap sizeMap){
		for(Bin b : bins) {
			if(b.size()>minSize) {b.calcContam(sizeMap);}
		}
		return toScoreString(toStats(bins, minSize), sizeMap.sum());
	}
	
	/**
	 * Generates a score string for BinStats objects.
	 * 
	 * @param bins List of BinStats to score
	 * @param totalSize Total size for calculating percentages
	 * @return Formatted score string
	 */
	private static String toScoreString(ArrayList<BinStats> bins, long totalSize){
		double compltScore=0, contamScore=0;
		double totalScore2=0;
		IntHashMap tidBins=new IntHashMap();
		for(BinStats bin : bins) {
			if(bin.taxid>0) {
				tidBins.increment(bin.taxid);
			}
			long contam=Math.round(bin.contam*bin.size);
			contamScore+=contam;
			compltScore+=Math.round(bin.complt*(bin.size-contam));
			double score=Math.max(0, bin.complt-5*bin.contam);
			totalScore2+=score*score;
		}
		String compS=String.format("%.3f", 100*compltScore/totalSize);
		String contamS=String.format("%.4f", 100*contamScore/totalSize);
		String totalS=String.format("%.2f", totalScore2);
		return "Complt:\t"+compS+"\tContam:\t"+contamS+"\tTotal:\t"+totalS;
	}
	
	/**
	 * Prints statistics about clean vs. contaminated bins.
	 * 
	 * @param bins List of BinStats to analyze
	 */
	public static void printCleanDirty(ArrayList<BinStats> bins) {
		long cleanBins=0, contamBins=0;
		long cleanContigs=0, contamContigs=0;
		long cleanSize=0, contamSize=0;
		long partialCleanSize=0, partialContamSize=0;
		long badContigs=0;
		for(BinStats bin : bins) {
			long contam=Math.round(bin.contam*bin.size);
			badContigs+=bin.badContigs;
			if(contam<1) {
				cleanBins++;
				cleanSize+=bin.size;
				cleanContigs+=bin.contigs;
			}else {
				contamBins++;
				contamSize+=bin.size;
				contamContigs+=bin.contigs;
				partialCleanSize+=(bin.size-contam);
				partialContamSize+=contam;
			}
		}
		outstream.println(QuickBin.formatString("Clean Bins", 29, cleanBins, contamBins));
		outstream.println(QuickBin.formatString("Dirty Bins", 29, contamBins, cleanBins));
		outstream.println(QuickBin.formatString("Clean Bin Bases", 29, cleanSize, contamSize));
		outstream.println(QuickBin.formatString("Dirty Bin Bases", 29, contamSize, cleanSize));
		outstream.println(QuickBin.formatString("Tainted Bases", 29, 
				partialCleanSize, cleanSize+contamSize-partialCleanSize));
		outstream.println(QuickBin.formatString("Contam Bases", 29, 
				partialContamSize, cleanSize+contamSize-partialContamSize));
		outstream.println("Bad Contigs:                 \t"+
				String.format("%.3f", badContigs*100.0/(cleanContigs+contamContigs)));
	}

	/**
	 * Loads bins using a single thread, typically for smaller datasets.
	 * 
	 * @param in List of input file paths to process
	 * @return List of BinStats objects representing the bins
	 */
	public static ArrayList<BinStats> loadST(ArrayList<String> in){
		ArrayList<BinStats> bins=new ArrayList<BinStats>(in.size());
		for(String s : in) {
			final BinStats bs;
			Cluster clust=loadCluster(s);
			calcContam(s, clust);
			bs=new BinStats(clust, ReadWrite.stripToCore(s));
			if(runQuickClade) {bs.taxid=callTax(clust);}
			if(callGenes) {
				callGenes(clust, GeneTools.gCaller, bs);
			}else if(gffMap!=null) {
				annotate(clust, gffMap, bs);
			}
			
			bins.add(bs);
		}
		return bins;
	}
	
	/**
	 * Loads bins using multiple threads for improved performance.
	 * 
	 * @param in List of input file paths to process
	 * @return List of BinStats objects representing the bins
	 */
	public ArrayList<BinStats> loadMT(ArrayList<String> in){
		//Do anything necessary prior to processing
		ArrayList<BinStats> bins=new ArrayList<BinStats>(in.size());
		
		//Determine how many threads may be used
		int threads=Shared.threads();
		if(threads>16) {threads=Tools.mid(16, threads/2, 32);}
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(in, bins, i, threads));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		this.success&=!success;
		
		//Do anything necessary after processing
		return bins;
	}

	/**
	 * Loads a cluster (bin) from a FASTA file.
	 * Reads contigs, calculates GC content, and processes taxonomy IDs.
	 * 
	 * @param fname Input file path
	 * @return Cluster object containing the bin's contigs
	 */
	static Cluster loadCluster(String fname) {
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			cris.start();
		}
		Cluster c=new Cluster(0);
		c.tetramers=new int[0];
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//ln!=null prevents a compiler potential null access warning
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					//  *********  Process reads here  *********
					Contig a=new Contig(r1.name(), r1.bases, (int)r1.numericID);
					for(byte b : a.bases) {
						int x=AminoAcid.baseToNumber[b];
						a.gcSum+=(x==1 || x==2) ? 1 : 0;
					}
					int tid=BinObject.parseTaxID(a.name);
					a.taxid=a.labelTaxid=tid;
					String key=ContigRenamer.toShortName(a.name);
					if(camiMap!=null) {
						Integer camiTid=camiMap.get(key);
						a.labelTaxid=(camiTid==null ? 0 : camiTid.intValue());
					}
					if(covMap!=null) {
						FloatList fl=covMap.get(key);
						if(fl!=null) {
							for(int i=0; i<fl.size; i++) {a.setDepth(fl.get(i), i);}
						}
					}
					c.add(a);
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		ReadWrite.closeStream(cris);
		return c;
	}
	
	/**
	 * Calculates contamination and completeness for a bin.
	 * Uses external quality assessment tools (CheckM, EukCC) if available.
	 * 
	 * @param fname Input file name for looking up external quality assessments
	 * @param c Cluster object to calculate metrics for
	 */
	static void calcContam(String fname, Cluster c) {
		fname=new File(fname).getName();
		String core=ReadWrite.stripToCore(fname);
		CCLine dummy=new CCLine(0, 0);
		CCLine checkm=(checkMMap==null ? null : checkMMap.get(core));
		CCLine eukcc=(eukCCMap==null ? null : eukCCMap.get(core));
		assert((checkMMap==null) == (checkm==null)) : checkm;
		if(checkm==null && eukcc==null) {
			c.calcContam(sizeMap);
			return;
		}
		if(checkm==null) {checkm=dummy;}
		if(eukcc==null) {eukcc=dummy;}
		CCLine best=(checkm.completeness>=eukcc.completeness ? checkm : eukcc);
		c.completeness=best.completeness;
		c.contam=best.contam;
	}
	
	/**
	 * Converts a collection of bins to BinStats objects for reporting.
	 * Optionally performs taxonomy assignment and gene calling.
	 * 
	 * @param bins Collection of bins to convert
	 * @param minSize Minimum size threshold for bins to include
	 * @return List of BinStats objects
	 */
	static ArrayList<BinStats> toStats(Collection<? extends Bin> bins, int minSize) {
		ArrayList<BinStats> list=new ArrayList<BinStats>();
		for(Bin b : bins) {
			if(b.size()>=minSize) {
				BinStats bs=new BinStats(b, b.name());
				if(runQuickClade) {bs.taxid=callTax(b);}
				if(callGenes) {
					callGenes(b, GeneTools.gCaller, bs);
				}else if(gffMap!=null) {
					annotate(b, gffMap, bs);
				}
				list.add(bs);
			}
		}
		return list;
	}
	
	/**
	 * Generates a detailed report file for all bins with quality metrics and annotations.
	 * 
	 * @param bins Collection of bins to report
	 * @param minSize Minimum size threshold for bins to include
	 * @param fname Output file path
	 */
	static void printClusterReport(Collection<? extends Bin> bins, int minSize, String fname) {
		ArrayList<BinStats> list=toStats(bins, minSize);
		printClusterReport(list, minSize, fname);
	}
	
	static void printClusterReport(ArrayList<BinStats> bins, int minSize, String fname) {
		if(fname==null) {return;}
		Collections.sort(bins);
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		String header="#Bin\tSize\tContigs\tGC\tDepth\tMinDepth\tMaxDepth";
		header+="\tCompleteness\tContam\tTaxID\tType";
		if(callGenes || gffFile!=null) {header+="\t16S\t18S\t23S\t5S\ttRNA\tCDS\tCDSLen";}
		if(BinObject.tree!=null) {header+="\tLineage";}
		bsw.println(header);
		for(BinStats b : bins) {
			if(b.size>=minSize) {
				bsw.printt(b.name).printt(b.size).printt(b.contigs);
				bsw.printt(b.gc, 3).printt(b.depth, 2);
				bsw.printt(b.minDepth, 2).printt(b.maxDepth, 2);
				bsw.printt(b.complt, 5).printt(b.contam, 5);
				bsw.printt(b.taxid).print(b.type(useRNA));
				
				if(callGenes || gffFile!=null) {
					bsw.tab().printt(b.r16Scount).printt(b.r18Scount);
					bsw.printt(b.r23Scount).printt(b.r5Scount);
					bsw.printt(b.trnaCount);
					bsw.printt(b.cdsCount).print(b.cdsLength);
				}
				
				if(BinObject.tree!=null) {
					bsw.tab().print(b.lineage!=null ? b.lineage : Clade.lineage(b.taxid));
				}
				bsw.println();
			}
		}
		bsw.poison();
	}
	
	/**
	 * Prints distribution of taxonomic classifications at different levels.
	 * 
	 * @param bins List of BinStats
	 * @param outstream Output stream to print to
	 */
	static void printTaxLevels(ArrayList<BinStats> bins, PrintStream outstream) {
		outstream.println("Unique Taxa Counts:");
		outstream.println("Level         \tTotal\tMQ\tHQ");
		for(int i=TaxTree.DOMAIN; i>=TaxTree.SPECIES; i--) {
			outstream.print(Tools.padRight(TaxTree.levelToString(i), 14));
			outstream.print("\t"+levelMaps[i].size());
			outstream.print("\t"+levelMapsMQ[i].size());
			outstream.print("\t"+levelMapsHQ[i].size());
			outstream.println();
		}
	}
	
	/**
	 * Prints quality distribution of bins (UHQ, VHQ, HQ, MQ, LQ, VLQ).
	 * 
	 * @param bins Collection of bins
	 * @param minSize Minimum size threshold for bins to include
	 * @param useRNA Whether to require rRNA genes for HQ designation
	 * @param outstream Output stream to print to
	 */
	static void printBinQuality(Collection<? extends Bin> bins, int minSize, boolean useRNA, 
			PrintStream outstream) {
		ArrayList<BinStats> list=toStats(bins, minSize);
		printBinQuality(list, minSize, useRNA, outstream);
	}
	
	static void printBinQuality(ArrayList<BinStats> bins, int minSize, boolean useRNA, 
			PrintStream outstream) {
		long uhq=0, uhqINC=0, uhqCON=0;
		long vhq=0, vhqINC=0, vhqCON=0;
		long hq=0, hqINC=0, hqCON=0;
		long mq=0, mqINC=0, mqCON=0;
		long lq=0, lqINC=0, lqCON=0;
		long vlq=0, vlqINC=0, vlqCON=0;

		long uhqSize=0;
		long vhqSize=0;
		long hqSize=0;
		long mqSize=0;
		long lqSize=0;
		long vlqSize=0;
		
		for(BinStats b : bins) {
			final long size=b.size;
			final float comp=b.complt, contam=b.contam;
			if(size>=minSize) {
				if(contam<=0.05f && comp>=0.9f && (!useRNA || (b.r16Scount>0 && b.r23Scount>0 && b.trnaCount>=18))) {
					hq++;
					hqSize+=size;
					if(comp>=0.99f && contam<=0.01f) {
						uhq++;
						uhqSize+=size;
						if(comp<1) {uhqINC++;}
						if(contam>0) {uhqCON++;}
					}else if(comp>=0.95f && contam<=0.02f) {
						vhq++;
						vhqSize+=size;
						if(comp<.99f) {vhqINC++;}
						if(contam>0.01f) {vhqCON++;}
					}else {
						if(comp<.95f) {hqINC++;}
						if(contam>0.02f) {hqCON++;}
					}
				}else if(contam<0.10f && comp>=0.5f) {
					mq++;
					mqSize+=size;
					if(comp<.90f) {mqINC++;}
					if(contam>0.05f) {mqCON++;}
				}else {
					lq++;
					lqSize+=size;
					if(contam>0.20f || comp<0.20f) {//vlq
						vlq++;
						vlqSize+=size;
						if(comp<0.2f) {vlqINC++;}
						if(contam>0.2f) {vlqCON++;}
					}else {//lq, not vlq
						if(comp<0.5f) {lqINC++;}
						if(contam>0.1f) {lqCON++;}
					}
				}
			}
		}
		//Make sets inclusive
		vhq+=uhq;
		vhqSize+=uhqSize;
		
		outstream.println("Quality\tBins\tIncomp\tContam\tBases");
		outstream.println("UHQ\t"+uhq+"\t"+uhqINC+"\t"+uhqCON+"\t"+uhqSize);
		outstream.println("VHQ\t"+vhq+"\t"+vhqINC+"\t"+vhqCON+"\t"+vhqSize);
		outstream.println("HQ\t"+hq+"\t"+hqINC+"\t"+hqCON+"\t"+hqSize);
		outstream.println("MQ\t"+mq+"\t"+mqINC+"\t"+mqCON+"\t"+mqSize);
		outstream.println("LQ\t"+lq+"\t"+lqINC+"\t"+lqCON+"\t"+lqSize);
		outstream.println("VLQ\t"+vlq+"\t"+vlqINC+"\t"+vlqCON+"\t"+vlqSize);
		String hqm=""+(hq+mq/4f);
		if(hqm.endsWith(".0")) {hqm=hqm.substring(0, hqm.length()-2);}
		outstream.println("HQ+MQ/4\t"+hqm+"\t\t\t"+(hqSize+mqSize/4));
	}
	
	/**
	 * Prints L90 statistics for a collection of bins.
	 * 
	 * @param bins Collection of bins
	 * @param basesLoaded Total size of the reference dataset
	 */
	static void printL90FromBins(Collection<? extends Bin> bins, long basesLoaded) {
		LongList sizes=new LongList(bins.size());
		for(Bin b : bins) {
			sizes.add(b.size());
		}
		GradeBins.printL90(sizes, basesLoaded);
	}
	
	/**
	 * Prints L90 statistics for bin sizes.
	 * 
	 * @param bins Collection of BinStats
	 * @param basesLoaded Total size of the reference dataset
	 */
	static void printL90(Collection<BinStats> bins, long basesLoaded) {
		LongList sizes=new LongList(bins.size());
		for(BinStats b : bins) {
			sizes.add(b.size);
		}
		GradeBins.printL90(sizes, basesLoaded);
	}
	
	/**
	 * Prints L90 statistics for a list of sizes.
	 * Reports L01, L10, L20, L50, and L90 values.
	 * 
	 * @param list List of sizes
	 * @param basesLoaded Total size of the reference dataset
	 */
	static void printL90(LongList list, long basesLoaded) {
		long c99=(long)(0.99f*basesLoaded);
		long c95=(long)(0.95f*basesLoaded);
		long c90=(long)(0.90f*basesLoaded);
		long c80=(long)(0.80f*basesLoaded);
		long c75=(long)(0.75f*basesLoaded);
		long c50=(long)(0.50f*basesLoaded);
		long c40=(long)(0.40f*basesLoaded);
		long c30=(long)(0.30f*basesLoaded);
		long c25=(long)(0.25f*basesLoaded);
		long c20=(long)(0.20f*basesLoaded);
		long c10=(long)(0.10f*basesLoaded);
		long c05=(long)(0.05f*basesLoaded);
		long c01=(long)(0.01f*basesLoaded);
		
		list.sort();
		list.reverse();
		long prev=0, sum2=0;
		for(int i=0; i<list.size(); i++) {
			long size=list.get(i);
			prev=sum2;
			sum2+=size;
			int num=i+1;

			if(sum2>=c01 && prev<c01) {System.err.println("L01: "+size+"\t"+"N01: "+num);}
//			if(sum2>=c05 && prev<c05) {System.err.println("L05: "+size+"\t"+"N05: "+num);}
			if(sum2>=c10 && prev<c10) {System.err.println("L10: "+size+"\t"+"N10: "+num);}
			if(sum2>=c20 && prev<c20) {System.err.println("L20: "+size+"\t"+"N20: "+num);}
//			if(sum2>=c25 && prev<c25) {System.err.println("L25: "+size+"\t"+"N25: "+num);}
//			if(sum2>=c30 && prev<c30) {System.err.println("L30: "+size+"\t"+"N30: "+num);}
//			if(sum2>=c40 && prev<c40) {System.err.println("L40: "+size+"\t"+"N40: "+num);}
			if(sum2>=c50 && prev<c50) {System.err.println("L50: "+size+"\t"+"N50: "+num);}
//			if(sum2>=c75 && prev<c75) {System.err.println("L75: "+size+"\t"+"N75: "+num);}
//			if(sum2>=c80 && prev<c80) {System.err.println("L80: "+size+"\t"+"N80: "+num);}
			if(sum2>=c90 && prev<c90) {System.err.println("L90: "+size+"\t"+"N90: "+num);}
//			if(sum2>=c95 && prev<c95) {System.err.println("L95: "+size+"\t"+"N95: "+num);}
//			if(sum2>=c99 && prev<c99) {System.err.println("L99: "+size+"\t"+"N99: "+num);}
		}
	}
	
	/**
	 * Creates a size map from a reference FASTA file.
	 * Maps taxonomy IDs to total sequence size.
	 * 
	 * @param fname Reference FASTA file path
	 * @return Map of taxonomy IDs to sequence sizes
	 */
	IntLongHashMap makeSizeMap(String fname) {
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		
		final ConcurrentReadInputStream cris;
		cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
		cris.start();
		
		IntLongHashMap map=new IntLongHashMap();
		countMap=new IntHashMap();
		long sizeSum=0, contigSum=0;
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//ln!=null prevents a compiler potential null access warning
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					readsProcessed++;
					basesProcessed+=r.length();
					sizeSum+=r.length();
					contigSum++;
					
					//  *********  Process reads here  *********
					int tid=BinObject.parseTaxID(r.id);
					long ret=map.increment(tid, r.length());
					countMap.increment(tid);
					if(ret==r.length() && tid>0) {taxIDsIn++;}
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		if(totalSize==0) {totalSize=sizeSum;}
		if(totalContigs==0) {totalContigs=contigSum;}
		ReadWrite.closeStream(cris);
		return map;
	}
	
	/**
	 * Loads taxonomy and size information from a tab-delimited file.
	 * 
	 * @param fname Input file path
	 * @return Map of taxonomy IDs to sequence sizes
	 */
	private IntLongHashMap loadTaxIn(String fname) {
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		LineParser1 lp=new LineParser1('\t');
		
		IntLongHashMap map=new IntLongHashMap();
		countMap=new IntHashMap();
		long sizeSum=0, contigSum=0;
		for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()) {
			for(byte[] line : ln) {
				lp.set(line);
				if(!lp.startsWith('#')) {
					int tid=lp.parseInt(0);
					long size=lp.parseLong(1);
					int contigs=lp.parseInt(2);
					long ret=map.increment(tid, size);
					countMap.increment(tid, contigs);
					sizeSum+=size;
					contigSum+=contigs;
					if(ret==size && tid>0) {taxIDsIn++;}
				}
			}
		}
		if(totalSize==0) {totalSize=sizeSum;}
		if(totalContigs==0) {totalContigs=contigSum;}
		bf.close();
		return map;
	}
	
	/**
	 * Writes taxonomy and size information to a tab-delimited file.
	 * 
	 * @param fname Output file path
	 * @param sizeMap Map of taxonomy IDs to sequence sizes
	 * @param countMap Map of taxonomy IDs to contig counts
	 */
	private void writeTaxOut(String fname, IntLongHashMap sizeMap, IntHashMap countMap) {
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, false);
		bsw.print("#taxID\tSize\tContigs\n");
		int[] tids=sizeMap.toArray();
		Arrays.sort(tids);
		for(int tid : tids) {
			bsw.print(tid).tab().print(sizeMap.get(tid)).tab().println(countMap.get(tid));
		}
		bsw.poison();
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads CheckM quality assessment results.
	 * 
	 * @param fname CheckM output file or directory
	 * @return Map of bin names to quality metrics
	 */
	public static HashMap<String, CCLine> loadCheckM(String fname){
		if(fname==null) {return null;}
		File f=new File(fname);
		if(f.isDirectory()) {
			if(!fname.endsWith("/")) {fname=fname+"/";}
			fname=fname+"quality_report.tsv";
		}
		ArrayList<byte[]> lines=ByteFile.toLines(fname);
		HashMap<String, CCLine> map=new HashMap<String, CCLine>(lines.size());
		LineParser1 lp=new LineParser1('\t');
		for(byte[] line : lines) {
			lp.set(line);
			if(!lp.startsWith("Name\t")) {
				String name=ReadWrite.stripToCore(lp.parseString(0));
				float comp=lp.parseFloat(1)/100;
				float contam=lp.parseFloat(2)/100;
				comp=Tools.mid(0, 1, comp);
				contam=Tools.mid(0, 1, contam);
				assert(comp>=0 && comp<=1) : new String(line);
//				assert(contam>=0 && contam<=1) : new String(line);
//				long size=//unavailable
				CCLine cc=new CCLine(comp, contam);
				map.put(name, cc);
			}
		}
		return map;
	}
	
	/**
	 * Loads EukCC quality assessment results.
	 * 
	 * @param fname EukCC output file or directory
	 * @return Map of bin names to quality metrics
	 */
	public static HashMap<String, CCLine> loadEukCC(String fname){
		if(fname==null) {return null;}
		File f=new File(fname);
		if(f.isDirectory()) {
			if(!fname.endsWith("/")) {fname=fname+"/";}
			fname=fname+"eukcc.csv";
		}
		ArrayList<byte[]> lines=ByteFile.toLines(fname);
		HashMap<String, CCLine> map=new HashMap<String, CCLine>(lines.size());
		LineParser1 lp=new LineParser1('\t');
		for(byte[] line : lines) {
			lp.set(line);
			if(!lp.startsWith("bin\tcompleteness")) {
				String name=ReadWrite.stripToCore(lp.parseString(0));
				float comp=lp.parseFloat(1)/100;
				float contam=lp.parseFloat(2)/100;
				comp=Tools.mid(0, 1, comp);
				contam=Tools.mid(0, 1, contam);
				assert(comp>=0 && comp<=1) : new String(line);
//				assert(contam>=0 && contam<=1) : new String(line);
//				long size=//unavailable
				CCLine cc=new CCLine(comp, contam);
				map.put(name, cc);
			}
		}
		return map;
	}
	
	/**
	 * Loads CAMI taxonomic assignments.
	 * 
	 * @param fname CAMI output file
	 * @return Map of contig names to taxonomy IDs
	 */
	public static HashMap<String, Integer> loadCami(String fname) {
		if(fname==null) {return null;}
		LineParser1 lp=new LineParser1('\t');
		ArrayList<byte[]> lines=ByteFile.toLines(fname);
		HashMap<String, Integer> map=new HashMap<String, Integer>();
		for(byte[] line : lines) {
			if(!Tools.startsWith(line, '@')){
				lp.set(line);
				String name=lp.parseString(0);
				int taxID=lp.parseInt(2);
				map.put(name, taxID);
			}
		}
		return map;
	}
	
	/**
	 * Loads GTDB-Tk taxonomic classifications from a directory.
	 * 
	 * @param fname GTDB-Tk output directory or file
	 * @return Map of bin names to taxonomic lineages
	 */
	public static HashMap<String, Lineage> loadGTDBDir(String fname) {
		if(fname==null) {return null;}
		HashMap<String, Lineage> map=new HashMap<String, Lineage>();
		File f=new File(fname);
		if(f.isDirectory()) {
			if(!fname.endsWith("/") && !fname.endsWith("\\")) {fname=fname+"/";}
			String bac=fname+"gtdbtk.bac120.summary.tsv";
			String ar=fname+"gtdbtk.ar53.summary.tsv";
			int loaded=0;
			if(loadGTDBFile(bac, map)) {loaded++;}
			if(loadGTDBFile(ar, map)) {loaded++;}
			assert(loaded>0) : "Could not find "+bac+" or "+ar;
		}else {
			loadGTDBFile(fname, map);
		}
		return map;
	}
	
	/**
	 * Loads a specific GTDB-Tk classification file.
	 * 
	 * @param fname GTDB-Tk output file
	 * @param map Map to populate with taxonomic lineages
	 * @return True if file was successfully loaded
	 */
	public static boolean loadGTDBFile(String fname, HashMap<String, Lineage> map) {
		if(fname==null || !new File(fname).canRead()) {return false;}
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		LineParser1 lptab=new LineParser1('\t');
		LineParser1 lpsemi=new LineParser1(';');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(line[0]=='u' && Tools.startsWith(line, "user_genome	classification")) {continue;}
			lptab.set(line);
			GTDBLine gline=new GTDBLine(lptab, lpsemi);
			if(map.containsKey(gline.name)) {continue;}//Only one taxa per bin
			map.put(gline.name, new Lineage(gline.classification));
		}
		return true;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------          Accumulator         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Accumulates results from worker threads.
	 * Part of the Accumulator interface.
	 * 
	 * @param t ProcessThread to accumulate results from
	 */
	@Override
	public void accumulate(ProcessThread t) {
		success=(success && t.success);
	}

	/**
	 * Returns a lock for thread synchronization.
	 * Part of the Accumulator interface.
	 * 
	 * @return ReadWriteLock for synchronization
	 */
	@Override
	public ReadWriteLock rwlock() {
		return null;
	}

	/**
	 * Indicates whether processing completed successfully.
	 * Part of the Accumulator interface.
	 * 
	 * @return True if processing was successful
	 */
	@Override
	public boolean success() {
		return success;
	}
	
	/**
	 * Predicts the taxonomy of a bin using k-mer composition.
	 * Requires the QuickClade index to be loaded.
	 * 
	 * @param b Bin to classify
	 * @return Predicted taxonomy ID or -1 if classification failed
	 */
	static int callTax(Bin b) {
		Clade clade=new Clade(-1, -1, b.name());
		for(Contig c : b) {
			clade.add(c.bases, null);
		}
		clade.finish();
		ArrayList<Comparison> list=cladeIndex.findBest(clade);
		Comparison best=(list==null ? null : list.get(0));
		return best==null || best.ref==null ? -1 : best.ref.taxID;
	}
	
	/**
	 * Identifies genes in a bin including 16S/18S/23S/5S rRNA, tRNA, and protein-coding genes.
	 * Requires gene caller to be initialized.
	 * 
	 * @param b Bin to annotate
	 * @param gcall GeneCaller instance for gene prediction
	 * @param bs BinStats object to update with gene counts
	 */
	static void callGenes(Bin b, GeneCaller gcall, BinStats bs) {
		ArrayList<Read> reads=new ArrayList<Read>(b.numContigs());
		for(Contig c : b) {
			reads.add(new Read(c.bases, null, c.name, c.id()));
		}
		ArrayList<Orf> orfs=gcall.callGenes(reads);
		for(Orf o : orfs) {
			if(o.is16S()) {bs.r16Scount++;}
			if(o.is18S()) {bs.r18Scount++;}
			if(o.is23S()) {bs.r23Scount++;}
			if(o.is5S()) {bs.r5Scount++;}
			if(o.isTRNA()) {bs.trnaCount++;}
			if(o.isCDS()) {
				bs.cdsCount++;
				bs.cdsLength+=o.length();
			}
		}
	}
	
	/**
	 * Annotates a bin with genes from pre-loaded GFF files.
	 * Updates gene counts in the BinStats object.
	 * 
	 * @param b Bin to annotate
	 * @param map HashMap mapping contig names to GFF lines
	 * @param bs BinStats object to update with gene counts
	 */
	static void annotate(Bin b, HashMap<String, ArrayList<GffLine>> map, BinStats bs) {
//		System.err.println("Annotating "+b.name());
		for(Contig c : b) {
			String name=c.name;
			ArrayList<GffLine> lines=map.get(name);
			if(lines==null) {lines=map.get(ContigRenamer.toShortName(name));}
			if(lines==null) {continue;}
			for(GffLine line : lines) {
				final int type=line.prokType();
				if(type==ProkObject.r16S) {bs.r16Scount++;}
				else if(type==ProkObject.r18S) {bs.r18Scount++;}
				else if(type==ProkObject.r23S) {bs.r23Scount++;}
				else if(type==ProkObject.r5S) {bs.r5Scount++;}
				else if(type==ProkObject.tRNA) {bs.trnaCount++;}
				else if(type==ProkObject.CDS) {
					bs.cdsCount++;
					bs.cdsLength+=line.length();
				}else {
					System.err.println("No match for "+line);
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier.
	 * 
	 * Thread for parallel processing of bin files.
	 * Each thread handles a subset of the input files.
	 */
	static class ProcessThread extends Thread {
		
		ProcessThread(ArrayList<String> fnames_, ArrayList<BinStats> bins_, 
				int tid_, int threads_){
			fnames=fnames_;
			bins=bins_;
			tid=tid_;
			threads=threads_;
			gCallerT=(callGenes ? GeneTools.makeGeneCaller() : null);
		}
		
		@Override
		public void run() {
			for(int i=tid; i<fnames.size(); i+=threads) {
				String fname=fnames.get(i);
				Cluster clust=loadCluster(fname);
				calcContam(fname, clust);
				BinStats bs=new BinStats(clust, ReadWrite.stripToCore(fname));
				if(runQuickClade) {bs.taxid=callTax(clust);}
				
				if(callGenes) {
					callGenes(clust, gCallerT, bs);
				}else if(gffMap!=null) {
					annotate(clust, gffMap, bs);
				}
				synchronized(bins) {
					bins.add(bs);
				}
			}
			success=true;
		}
		
		private final ArrayList<String> fnames;
		private final ArrayList<BinStats> bins;
		private final int tid;
		private final int threads;
		private final GeneCaller gCallerT;
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Container for completeness and contamination values.
	 * Used for storing quality metrics from external tools.
	 */
	private static class CCLine {
		
		CCLine(float completeness_, float contam_) {
			this(completeness_, contam_, -1);
		}
		
		CCLine(float completeness_, float contam_, long size_) {
			completeness=completeness_;
			contam=contam_;
			size=size_;
			assert(completeness>=0 && completeness<=1) : completeness;
			assert(contam>=0 && contam<=1) : contam;
			assert(size>0 || size==-1);
		}
		
		public String toString() {return size+", "+completeness+", "+contam;}
		
		long size=-1;
		float completeness=-1;
		float contam=-1;
		
	}
	
	/*--------------------------------------------------------------*/
	
	// Input/output files
	private ArrayList<String> in;         // Input bin files
	private String taxIn;                 // Input taxonomy file
	private String taxOut;                // Output taxonomy file
	private String tax;                   // Combined taxonomy file (input or output)
	private String ref;                   // Reference genome file
	private String hist;                  // Output histogram file
	private String contamHist;            // Output contamination histogram file
	private String ccplot;                // Output completeness-contamination plot file
	private String checkMFile;            // CheckM output file
	private String eukCCFile;             // EukCC output file
	private String camiFile;              // CAMI classification file
	private String gtdbFile;              // GTDB-Tk classification file
	private static String cov;            // Coverage information file
	private static String gffFile;        // GFF annotation file
	private static String imgMapFile;     // IMG mapping file
	private static String spectraFile;    // QuickClade reference file

	// Output report file
	private String report;                // Detailed report output file

	// Data structures
	private LongList sizes;               // List of bin sizes
	private ArrayList<BinStats> bins;     // List of processed bins
	private double contamScore;           // Total contamination score
	private double compltScore;           // Total completeness score
	private int minSize;                  // Minimum bin size to include
	
	/*--------------------------------------------------------------*/

	// Processing flags
	private boolean loadMT;               // Use multi-threading for loading
	private static boolean runQuickClade; // Perform k-mer-based taxonomy assignment
	private static boolean useTree;       // Use taxonomy tree for classification
	private static boolean callGenes;     // Identify genes in bins
	private static boolean useRNA;        // Require RNA genes for high-quality designation

	// Maps
	private static IntLongHashMap sizeMap;    // Maps taxonomy IDs to genome sizes
	private static IntHashMap countMap;       // Maps taxonomy IDs to contig counts
	private static HashMap<String, CCLine> checkMMap;      // CheckM quality results
	private static HashMap<String, CCLine> eukCCMap;       // EukCC quality results
	private static HashMap<String, Integer> camiMap;       // CAMI taxonomy assignments
	private static HashMap<String, Lineage> gtdbMap;       // GTDB-Tk taxonomic lineages
	private static HashMap<String, ArrayList<GffLine>> gffMap;  // GFF annotations
	private static HashMap<String, FloatList> covMap;      // Coverage information

	// Taxonomy tracking
	private static IntHashMap[] levelMaps;    // Counts organisms at each taxonomic level
	private static IntHashMap[] levelMapsHQ;  // Counts high-quality organisms at each level
	private static IntHashMap[] levelMapsMQ;  // Counts medium-quality organisms at each level
	private static CladeIndex cladeIndex;     // Index for k-mer-based classification

	// Statistics
	private static long maxReads;         // Maximum reads to process
	private long readsProcessed;          // Reads processed so far
	private long basesProcessed;          // Bases processed so far
	private long totalSize;               // Total size of reference
	private long totalContigs;            // Total contigs in reference
	private long taxIDsIn;                // Unique taxonomy IDs in reference

	// Control flags
	boolean overwrite;                    // Overwrite existing output files
	boolean success;                      // Processing completed successfully
	
	/*--------------------------------------------------------------*/

	// Output control
	private static java.io.PrintStream outstream;  // Output stream
	public static boolean verbose;        // Verbose output flag
	
	/*--------------------------------------------------------------*/
	
}
