package prok;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.PrintStream;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;


import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.CompareGff;
import gff.GffLine;
import jgi.BBMerge;
import json.JsonObject;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;
import tracker.EntropyTracker;
import ml.CellNet;
import ml.CellNetParser;
import prok.CallGenesHelper.ContigStats;
import prok.CallGenesHelper.GeneQuad;
import prok.CallGenesHelper.TrueGeneData;

/**
 * This is the executable class for gene-calling.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class CallGenes extends ProkObject {
	
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
		CallGenes x=new CallGenes(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CallGenes(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, (args.length>40 ? null : getClass()), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			//outGff=parser.out1; // Possibly corrupting outGff with null
			maxReads=parser.maxReads;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program
		
		ffoutGff=FileFormat.testOutput(outGff, FileFormat.GFF, null, true, overwrite, append, ordered);
		ffoutAmino=FileFormat.testOutput(outAmino, FileFormat.FA, null, true, overwrite, append, ordered);
		ffout16S=FileFormat.testOutput(out16S, FileFormat.FA, null, true, overwrite, append, ordered);
		ffout18S=FileFormat.testOutput(out18S, FileFormat.FA, null, true, overwrite, append, ordered);
		ffnet = FileFormat.testOutput(netFile, FileFormat.BBNET, null, true, true, false, false); // Brandon: netFile is for CellNetParser
		
		if(ffoutGff!=null){
			assert(!ffoutGff.isSequence()) : "\nout is for gff files.  To output sequence, please use outa.";
		}
		if(ffoutAmino!=null){
			assert(!ffoutAmino.gff()) : "\nouta is for sequence data.  To output gff, please use out.";
		}
		if(ffout16S!=null){
			assert(!ffout16S.gff()) : "\nout16S is for sequence data.  To output gff, please use out.";
		}
		if(ffout18S!=null){
			assert(!ffout18S.gff()) : "\nout18S is for sequence data.  To output gff, please use out.";
		}
		
		if(geneHistFile==null){geneHistBins=0;}
		else{
			assert(geneHistBins>1) : "geneHistBins="+geneHistBins+"; should be >1";
			assert(geneHistDiv>=1) : "geneHistDiv="+geneHistDiv+"; should be >=1";
		}
		geneHist=geneHistBins>1 ? new long[geneHistBins] : null;

		if(netFile != null){ // Brandon: Load the network file if specified
			net0 = CellNetParser.load(netFile);
			assert(net0 != null) : netFile;
		} else {
			net0 = null;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

//			outstream.println(arg+", "+a+", "+b);
			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || 
					a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(b==null && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else if(a.equals("pgm") || a.equals("gm") || a.equals("model")){
				assert(b!=null);
				if(b.equalsIgnoreCase("auto") || b.equalsIgnoreCase("default")){
					b=Data.findPath("?model.pgm");
					pgmList.add(b);
				}else{
					Tools.addFiles(b, pgmList);
				}
			}else if(b==null && new File(arg).exists() && FileFormat.isPgmFile(arg)){
				pgmList.add(b);
			}else if(a.equals("outamino") || a.equals("aminoout") || a.equals("outa") || a.equals("outaa") || a.equals("aaout") || a.equals("amino")){
				outAmino=b;
			}else if(a.equalsIgnoreCase("out16s") || a.equalsIgnoreCase("16sout")){
				out16S=b;
			}else if(a.equalsIgnoreCase("out18s") || a.equalsIgnoreCase("18sout")){
				out18S=b;
			}
			else if(a.equals("truegenes")){  // Brandon: truegenes is the file containing true gene annotations
				trueGenesFile = b;
				try {
					// Call the new helper method and unpack the results
					TrueGeneData data = CallGenesHelper.loadTrueGeneSet(trueGenesFile);
					this.trueGeneSet = data.set();
					this.totalGffRows = data.totalRows();
					this.totalGffGeneRows = data.geneRows();
				} catch (IOException e) {
					System.err.println("CRITICAL ERROR: Failed to load true genes file: " + trueGenesFile);
					e.printStackTrace();
				}
			}
			else if(arg.equalsIgnoreCase("nofilter")){ // Brandon: arg for no filtering of candidates
    			nofilter = true;
			}else if(a.equals("net")){ // Brandon: arg for network input
				netFile = b;
			}else if (arg.equalsIgnoreCase("seq")) { // Brandon: seq is to output the sequence of the gene, rather than one-hot
				seq = true;
			}
			else if(a.equals("outgff") || a.equals("gffout") || a.equals("outg") || a.equals("gff")){
				outGff=b;
			}
			else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				//ReadWrite.verbose=verbose;
				GeneCaller.verbose=verbose;
			}else if(a.equalsIgnoreCase("ingff") || a.equalsIgnoreCase("gffin")){
				Tools.addFiles(b, inGffList);
			}
			
			else if(a.equals("json_out") || a.equalsIgnoreCase("json")){
				json_out=Parse.parseBoolean(b);
			}else if(a.equals("stats") || a.equalsIgnoreCase("outstats")){
				outStats=b;
			}else if(a.equals("hist") || a.equalsIgnoreCase("outhist") || a.equalsIgnoreCase("lengthhist") || a.equalsIgnoreCase("lhist") || a.equalsIgnoreCase("genehist")){
				geneHistFile=b;
			}else if(a.equals("bins")){
				geneHistBins=Integer.parseInt(b);
			}else if(a.equals("binlen") || a.equals("binlength") || a.equalsIgnoreCase("histdiv")){
				geneHistDiv=Integer.parseInt(b);
			}else if(a.equals("printzero") || a.equals("pz")){
				printZeroCountHist=Parse.parseBoolean(b);
			}
			
			else if(a.equals("merge")){
				merge=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("2pass") || a.equals("twopass")){
				if(Parse.parseBoolean(b)) {passes=2;}
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}

			else if(a.equals("extended") || a.equals("extendedstats") || a.equals("verbosestats")){
				extendedStats=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("setbias16s")) {
				GeneCaller.biases[r16S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbias18s")) {
				GeneCaller.biases[r18S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbias23s")) {
				GeneCaller.biases[r23S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbias5s")) {
				GeneCaller.biases[r5S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbiastRNA")) {
				GeneCaller.biases[tRNA]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbiasCDS")) {
				GeneCaller.biases[CDS]=Float.parseFloat(b);
			}
			
			else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}
			
			else if(a.equals("translate")){
				mode=TRANSLATE;
			}else if(a.equals("retranslate") || a.equals("detranslate")){
				mode=RETRANSLATE;
			}else if(a.equals("recode")){
				mode=RECODE;
			}
			
			else if(a.equalsIgnoreCase("minlen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("maxoverlapss") || a.equals("overlapss") || a.equals("overlapsamestrand") || a.equals("moss") || a.equalsIgnoreCase("maxOverlapSameStrand")){
				maxOverlapSameStrand=Integer.parseInt(b);
			}else if(a.equals("maxoverlapos") || a.equals("overlapos") || a.equals("overlapoppositestrand") || a.equals("moos") || a.equalsIgnoreCase("maxOverlapOppositeStrand")){
				maxOverlapOppositeStrand=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minStartScore")){
				minStartScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minStopScore")){
				minStopScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minInnerScore") || a.equalsIgnoreCase("minKmerScore")){
				minKmerScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minOrfScore") || a.equalsIgnoreCase("minScore")){
				minOrfScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minAvgScore")){
				minAvgScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("breakLimit")){
				GeneCaller.breakLimit=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("e1")){
				Orf.e1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("e2")){
				Orf.e2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("e3")){
				Orf.e3=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("f1")){
				Orf.f1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("f2")){
				Orf.f2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("f3")){
				Orf.f3=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("p0")){
				GeneCaller.p0=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p1")){
				GeneCaller.p1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p2")){
				GeneCaller.p2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p3")){
				GeneCaller.p3=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p4")){
				GeneCaller.p4=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p5")){
				GeneCaller.p5=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p6")){
				GeneCaller.p6=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("q1")){
				GeneCaller.q1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q2")){
				GeneCaller.q2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q3")){
				GeneCaller.q3=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q4")){
				GeneCaller.q4=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q5")){
				GeneCaller.q5=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("lookback")){
				GeneCaller.lookbackPlus=GeneCaller.lookbackMinus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("lookbackplus")){
				GeneCaller.lookbackPlus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("lookbackminus")){
				GeneCaller.lookbackMinus=Integer.parseInt(b);
			}
			
			else if(a.equalsIgnoreCase("compareto")){
				compareToGff=b;
			}
			
			else if(ProkObject.parse(arg, a, b)){}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		if(pgmList.isEmpty()){
			String b=Data.findPath("?model.pgm");
			pgmList.add(b);
		}
		for(int i=0; i<pgmList.size(); i++){
			String s=pgmList.get(i);
			if(s.equalsIgnoreCase("auto") || s.equalsIgnoreCase("default")){
				pgmList.set(i, Data.findPath("?model.pgm"));
			}
		}
		
		if(Shared.threads()<2){ordered=false;}
		assert(!fnaList.isEmpty()) : "At least 1 fasta file is required.";
		assert(!pgmList.isEmpty()) : "At least 1 pgm file is required.";
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		pgmList=Tools.fixExtension(pgmList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outGff, outAmino, out16S, out18S, outStats, geneHistFile)){
			outstream.println((outGff==null)+", "+outGff);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "
					+outGff+", "+outAmino+", "+out16S+", "+out18S+", "+outStats+", "+geneHistFile+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(pgmList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(outGff);
		foo.add(outAmino);
		foo.add(out16S);
		foo.add(out18S);
		foo.add(outStats);
		foo.add(geneHistFile);
		if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create read streams and process all data */
	void process(Timer t){


		try {
        // We assume the first FNA file is the reference for all contigs.
        // This is a typical use case.
			if (fnaList != null && !fnaList.isEmpty()) {
				String fastaFile = fnaList.get(0);
				System.err.println("Reading FASTA file: " + fastaFile);
				Map<String, String> contigSequences = CallGenesHelper.readFastaFile(fastaFile);
				this.contigMetrics = CallGenesHelper.calculateContigMetrics(contigSequences);
			}
			} catch (IOException e) {
				System.err.println("Error: Failed to read or process the FASTA file for contig metrics.");
				e.printStackTrace();
				// We can choose to exit or continue without the contig-level stats.
				// For now, we'll just print the error and continue.
				this.contigMetrics = new HashMap<>(); // Ensure it's not null
		}
			// Debugging output for contig IDs
		    /*System.err.println("--- Contig IDs found in FASTA map ---");
			if(this.contigMetrics != null) {
				for (String key : this.contigMetrics.keySet()) {
					System.err.println("MAP KEY: '" + key + "'");
				}
			}
			System.err.println("------------------------------------");
			*/
		final GeneModel pgm0=PGMTools.loadAndMerge(pgmList);
		
		if(call16S || call18S || call23S || calltRNA || call5S){
			loadLongKmers();
			loadConsensusSequenceFromFile(false, false);
		}
		
		ByteStreamWriter bsw=makeBSW(ffoutGff);
		if(bsw!=null){
			bsw.forcePrint("##gff-version 3\n");
		}

		ConcurrentReadOutputStream rosAmino=makeCros(ffoutAmino);
		ConcurrentReadOutputStream ros16S=makeCros(ffout16S);
		ConcurrentReadOutputStream ros18S=makeCros(ffout18S);

	
			
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsIn=genesOut=0;
		basesIn=bytesOut=0;
		
		for(int fnum=0; fnum<fnaList.size(); fnum++){
			final String fna=fnaList.get(fnum);
			String gffIn=(inGffList!=null && !inGffList.isEmpty()) ? inGffList.set(fnum, null) : null;
			//Create a read input stream
			final GeneModel pgm=makeMultipassModel(pgm0, fna, gffIn, passes/*, maxReads*/);
			
			final ConcurrentReadInputStream cris=makeCris(fna);
			
			//Process the reads in separate threads
			spawnThreads(cris, bsw, rosAmino, ros16S, ros18S, pgm);
			
			//Close the input stream
			errorState|=ReadWrite.closeStream(cris);
		}
		
		//Close the input stream
		errorState|=ReadWrite.closeStreams(null, rosAmino, ros16S, ros18S);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the output stream
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		outstream.println(Tools.linesBytesOut(readsIn, basesIn, genesOut, bytesOut, 8, false));
		outstream.println();
		
		if(json_out){
			printStatsJson(outStats);
		}else{
			printStats(outStats);
		}
		if(geneHistFile!=null){
			printHist(geneHistFile);
		}
		if(trueGenesFile!=null){
		outstream.println(); // Spacing for readability
		outstream.println("True Genes File: " + trueGenesFile);
		outstream.println("Total Rows detected in .gff = " + totalGffRows);
		outstream.println(".gff Gene Rows = " + totalGffGeneRows);
		outstream.println("Total Matches = " + totalMatches);
		}

		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		
		if(compareToGff!=null){
			if(compareToGff.equals("auto")){
				compareToGff=fnaList.get(0).replace(".fna", ".gff");
				compareToGff=compareToGff.replace(".fa", ".gff");
				compareToGff=compareToGff.replace(".fasta", ".gff");
			}
			CompareGff.main(new String[] {outGff, compareToGff});
		}
	}
	
	private void printStats(String fname){
		if(fname==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, false);
		bsw.start();
		
		if(callCDS){
			if(extendedStats || verbose) {
				bsw.println("Gene Starts Made:     \t "+Tools.padLeft(geneStartsMade, 12));
				bsw.println("Gene Stops Made:      \t "+Tools.padLeft(geneStopsMade, 12));
				bsw.println("Gene Starts Retained: \t "+Tools.padLeft(geneStartsRetained, 12));
				bsw.println("Gene Stops Retained:  \t "+Tools.padLeft(geneStopsRetained, 12));
			}
			bsw.print("Approx Coding Fraction:\t ");
			bsw.print(Tools.padLeft(Tools.format("%.4f", Tools.min(1.0, stCdsPass.lengthSum/(double)basesIn)), 12)).nl();
			bsw.println("CDS Out:              \t "+Tools.padLeft(geneStartsOut, 12));
		}
		if(call16S){bsw.println("16S Out:              \t "+Tools.padLeft(r16SOut, 12));}
		if(call18S){bsw.println("18S Out:              \t "+Tools.padLeft(r18SOut, 12));}
		if(call23S){bsw.println("23S Out:              \t "+Tools.padLeft(r23SOut, 12));}
		if(call5S){bsw.println("5S Out:               \t "+Tools.padLeft(r5SOut, 12));}
		if(calltRNA){bsw.println("tRNA Out:             \t "+Tools.padLeft(tRNAOut, 12));}
		
		if(extendedStats || verbose) {
			if(callCDS){
				bsw.println();
				bsw.println("All ORF Stats:");
				bsw.print(stCds.toString());

				//			bsw.println();
				//			bsw.println("Retained ORF Stats:");
				//			bsw.print(stCds2.toString());

				bsw.println();
				bsw.println("Called ORF Stats:");
				stCdsPass.genomeSize=basesIn;
				bsw.print(stCdsPass.toString());
			}

			if(call16S && r16SOut>0){
				bsw.println();
				bsw.println("Called 16S Stats:");
				bsw.print(st16s.toString());
			}
			if(call23S && r23SOut>0){
				bsw.println();
				bsw.println("Called 23S Stats:");
				bsw.print(st23s.toString());
			}
			if(call5S && r5SOut>0){
				bsw.println();
				bsw.println("Called 5S Stats:");
				bsw.print(st5s.toString());
			}
			if(call18S && r18SOut>0){
				bsw.println();
				bsw.println("Called 18S Stats:");
				bsw.print(st18s.toString());
			}
			if(calltRNA && tRNAOut>0){
				bsw.println();
				bsw.println("Called tRNA Stats:");
				bsw.print(sttRNA.toString());
			}
		}
		bsw.poisonAndWait();
	}
	
	private void printStatsJson(String fname){
		if(fname==null){return;}
		
		JsonObject outer=new JsonObject();
		
		{
			JsonObject jo=new JsonObject();
			if(callCDS){
				jo.add("Gene Starts Made", geneStartsMade);
				jo.add("Gene Stops Made", geneStopsMade);
				jo.add("Gene Starts Retained", geneStartsRetained);
				jo.add("Gene Stops Retained", geneStopsRetained);
				jo.add("CDS Out", geneStartsOut);
			}
			if(call16S){jo.add("16S Out", r16SOut);}
			if(call18S){jo.add("18S Out", r18SOut);}
			if(call23S){jo.add("23S Out", r23SOut);}
			if(call5S){jo.add("5S Out", r5SOut);}
			if(calltRNA){jo.add("tRNA Out", tRNAOut);}
			outer.add("Overall", jo);
		}
		
		if(callCDS){
			outer.add("All ORF Stats", stCds.toJson());
			outer.add("Retained ORF Stats", stCds2.toJson());
			stCdsPass.genomeSize=basesIn;
			outer.add("Called ORF Stats", stCdsPass.toJson());
		}

		if(call16S){outer.add("Called 16S Stats", st16s.toJson());}
		if(call18S){outer.add("Called 18S Stats", st18s.toJson());}
		if(call23S){outer.add("Called 23S Stats", st23s.toJson());}
		if(call5S){outer.add("Called 5S Stats", st5s.toJson());}
		if(calltRNA){outer.add("Called tRNA Stats", sttRNA.toJson());}
		

		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, false);
		bsw.start();
		bsw.println(outer.toText());
		bsw.poisonAndWait();
	}
	
	private void printHist(String fname){
		if(fname==null || geneHist==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, false);
		bsw.start();
		long sum=shared.Vector.sum(geneHist);
		double mean=Tools.averageHistogram(geneHist)*geneHistDiv;
		int median=Tools.medianHistogram(geneHist)*geneHistDiv;
		double std=Tools.standardDeviationHistogram(geneHist)*geneHistDiv;
		bsw.println("#Gene Length Histogram");
		bsw.print("#Genes:\t").println(sum);
		bsw.print("#Mean:\t").println(mean, 4);
		bsw.print("#Median:\t").println(median);
		bsw.print("#STDDev:\t").println(std, 4);
		bsw.print("#Length\tCount\n");
		long cum=0;
		for(int i=0; i<geneHist.length && cum<sum; i++){
			int len=i*geneHistDiv;
			long count=geneHist[i];
			cum+=count;
			if(count>0 || printZeroCountHist){
				bsw.print(len).tab().println(count);
			}
		}
		bsw.poisonAndWait();
	}
	
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ByteStreamWriter bsw, 
			ConcurrentReadOutputStream rosAmino, ConcurrentReadOutputStream ros16S, ConcurrentReadOutputStream ros18S, GeneModel pgm){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, rosAmino, ros16S, ros18S, pgm, minLen, i, this.contigMetrics,this.nofilter));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		
		// After all threads join, aggregate and print the total foundCount
		printTotalMatchCount(alpt);
	}

	/**
	 * Aggregates and prints the total number of found genes from all threads.
	 */
	private void printTotalMatchCount(ArrayList<ProcessThread> alpt) {
		int totalMatch = 0;
		for (ProcessThread pt : alpt) {
			totalMatch += pt.getMatchCount();
		}
		totalMatches = totalMatch;
	}
	
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
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			genesOut+=pt.genesOutT;
			bytesOut+=pt.bytesOutT;
			
			geneStopsMade+=pt.caller.geneStopsMade;
			geneStartsMade+=pt.caller.geneStartsMade;
			geneStartsRetained+=pt.caller.geneStartsRetained;
			geneStopsRetained+=pt.caller.geneStopsRetained;
			geneStartsOut+=pt.caller.geneStartsOut;

			r16SOut+=pt.caller.r16SOut;
			r18SOut+=pt.caller.r18SOut;
			r23SOut+=pt.caller.r23SOut;
			r5SOut+=pt.caller.r5SOut;
			tRNAOut+=pt.caller.tRNAOut;
			
			stCds.add(pt.caller.stCds);
			stCds2.add(pt.caller.stCds2);
//			stCdsPass.add(pt.caller.stCdsPass);
			
			for(int i=0; i<trackers.length; i++){
				trackers[i].add(pt.caller.trackers[i]);
			}
			
			if(geneHist!=null){Tools.add(geneHist, pt.geneHistT);}
			
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Multipass           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static GeneModel makeMultipassModel(GeneModel pgm0, String fna, String gff, int passes/*, long maxReads*/) {
		if(passes<2) {return pgm0;}
		ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1/*maxReads*/);
		return makeMultipassModel(pgm0, reads, gff, passes);
	}
	
	public static GeneModel makeMultipassModel(GeneModel pgm0, ArrayList<Read> reads, String gff, int passes) {
		GeneModel pgm=pgm0;
		for(int i=1; i<passes; i++) {
			pgm=runOnePass(reads, gff, pgm);
			if(i==1 && passes==2) {//only pass for 2-pass
				pgm=PGMTools.mix(0.50, 0.10, 0.10, true, pgm0, pgm);
			}else if(i==passes-1) {//final pass for 3+ passes
				pgm=PGMTools.mix(0.50, 0.12, 0.12, true, pgm0, pgm);
			}else if(i==1 && passes>2) {//first pass for 3+ passes
				pgm=PGMTools.mix(0.10, 0.01, 0.01, true, pgm0, pgm);
			}else {//middle pass for 4+ passes
				pgm=PGMTools.mix(0.20, 0.02, 0.02, true, pgm0, pgm);
			}
		}
		return pgm;
	}
	
	/** This needs a pgm OR a gff, not both */
	public static GeneModel runOnePass(String fna, String gff, GeneModel pgm0/*, long maxReads*/) {//TODO: Make this multithreaded
		ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1/*maxReads*/);
		return runOnePass(reads, gff, pgm0);
	}
	
	/** This needs a pgm OR a gff, not both */
	public static GeneModel runOnePass(ArrayList<Read> reads, String gff, GeneModel pgm0) {//TODO: Make this multithreaded
		GeneCaller caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, 
				minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm0);
		
		final ArrayList<GffLine> cds, rrna, trna;
		if(gff!=null && !"null".equalsIgnoreCase(gff)) {
			ArrayList<GffLine>[] allGffLines=GffLine.loadGffFileByType(gff, "CDS,rRNA,tRNA", true);
			cds=allGffLines[0];
			rrna=allGffLines[1];
			trna=allGffLines[2];
//			System.err.println("Loaded "+cds.size()+" CDSs");
		}else {
			ArrayList<Orf> orfs=caller.callGenes(reads);
			ArrayList<GffLine>[] allGffLines=Orf.toGffLinesByType(orfs, "CDS,rRNA,tRNA");
			cds=allGffLines[0];
			rrna=allGffLines[1];
			trna=allGffLines[2];
//			System.err.println("Called "+cds.size()+" CDSs");
		}
		
		GeneModel pgm=new GeneModel(true);
		pgm.process(reads, cds, rrna, trna);
		for(StatsContainer sc : pgm.allContainers) {
			sc.calculate();//Not sure if this is needed...
		}
		return pgm;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	private ConcurrentReadOutputStream makeCros(FileFormat ff){
		if(ff==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(4, 64, (Shared.threads()*2)/3) : 4);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	// Utility method to open a BufferedReader for plain or gzipped files
	private static BufferedReader openBufferedReader(String filename) throws IOException {
	    if (filename.endsWith(".gz")) {
	        return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename))));
	    } else {
	        return new BufferedReader(new FileReader(filename));
	    }
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {

        final Map<String, ContigStats> contigMetrics;
        final boolean nofilter;
        // FIXED: Declare thread-local network and feature vector fields
        private CellNet net;
        private float[] vec;
        
        //Constructor
        ProcessThread(final ConcurrentReadInputStream cris_, final ByteStreamWriter bsw_, 
                ConcurrentReadOutputStream rosAmino_, ConcurrentReadOutputStream ros16S_, ConcurrentReadOutputStream ros18S_, 
                GeneModel pgm_, final int minLen, final int tid_,
                final Map<String, ContigStats> contigMetrics_,
                final boolean nofilter_){
            cris=cris_;
            bsw=bsw_;
            rosAmino=rosAmino_;
            ros16S=ros16S_;
            ros18S=ros18S_;
            pgm=pgm_;
            tid=tid_;
            geneHistT=(geneHistBins>1 ? new long[geneHistBins] : null);
            caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, 
                    minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm);
            this.entropyTracker = new EntropyTracker(5,50,false);
            this.contigMetrics = contigMetrics_;
            this.nofilter = nofilter_; 
        }
        
        //Called by start()
        @Override
        public void run(){
            // FIXED: This block initializes the network and vector for each thread.
            // It MUST be present and run before processInner().
            if (net0 != null) {
                net = net0.copy(false);
                vec = new float[net.numInputs()];
            }
            
            //Process the reads
            processInner();
            
            //Indicate successful exit status
            success=true;
        }
        
        /** Iterate through the reads */
        void processInner(){
            
            //Grab the first ListNum of reads
            ListNum<Read> ln=cris.nextList();

            //As long as there is a nonempty read list...
            while(ln!=null && ln.size()>0){
                processList(ln);
                //Fetch a new list
                ln=cris.nextList();
            }

            //Notify the input stream that the final list was used
            if(ln!=null){
                cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
            }
        }
		
		// Brandon: 6/24/25
		/** Updated processList to handle CallGenesHelper.java */
		// 7/23/25: Updated to handle neural network scoring and filtering.
		// In CallGenes.java, inside the ProcessThread inner class
        void processList(ListNum<Read> ln) {
            final ArrayList<Read> reads = ln.list;

            // Loop through each contig (Read) in the list
            for (int idx = 0; idx < reads.size(); idx++) {
                final Read r1 = reads.get(idx);
                
                readsInT += r1.pairCount();
                basesInT += r1.length() + r1.mateLength();
                
                // Step 1: Get all potential gene candidates from the GeneCaller.
                // The 'nofilter' flag is handled inside processRead -> caller.callGenes.
                ArrayList<Orf> orfList = processRead(r1);
                ArrayList<Orf> passedOrfs;

                // Step 2: Filter candidates with the neural network if it's provided.
                if (net != null && !this.nofilter) {
                    // NN INFERENCE MODE: Score ORF candidates and keep only those that pass the cutoff.
                    passedOrfs = new ArrayList<>();
                    for (Orf orf : orfList) {
                        // Non-CDS features (tRNA, rRNA) are not scored by the network, so they pass automatically.
                        if (orf.type != CDS) {
                            passedOrfs.add(orf);
                            continue;
                        }
                        
                        // Generate feature vector, score it, and check if it passes the filter.
                        float score = scoreOrf(orf, r1);
                        boolean pass = (score >= cutoff) == highpass;
                        
                        if (pass) {
                            orf.orfScore = score; // Optionally update the ORF's score to the NN score.
                            passedOrfs.add(orf);
                        }
                    }
                } else {
                    // NO-NETWORK MODE: Keep all candidates. This branch is also used for 'nofilter'.
                    passedOrfs = orfList;
                }
                
                // Step 3: Process the final list of ORFs.
                if (passedOrfs == null || passedOrfs.isEmpty()) {
                    continue;
                }
                
                genesOutT += passedOrfs.size();

                // This buffer will hold all the GFF lines for the current contig.
                final ByteBuilder bb = new ByteBuilder();

                // Loop through the final set of ORFs (either filtered or unfiltered).
                for (Orf orf : passedOrfs) {
                    if (bsw != null) {
                        orf.appendGff(bb); // Generate the standard GFF line.
                        
                        // If not using the network, append the feature vector for training data generation.
                        if (net == null && orf.type == CDS) {
                            String vectorValues = CallGenesHelper.generateFeatureVector(
                                orf, r1, this.contigMetrics.get(r1.id), this.entropyTracker,
                                CallGenes.this.trueGeneSet, CallGenes.this.seq, false
                            );
                            
                            if (vectorValues.endsWith("\t1")) {
                                this.matchCount++;
                            }
                            
                            // Trim the newline from appendGff before adding the vector.
                            if (bb.length() > 0 && bb.get(bb.length() - 1) == '\n') {
                                bb.trimLast(1);
                            }
                            
                            // Append the vector attribute.
                            bb.append(";VECTOR=");
                            for (int i = 0; i < vectorValues.length(); i++) {
                                char c = vectorValues.charAt(i);
                                bb.append(c == '\t' ? ',' : c);
                            }
                            bb.nl();
                        }
                    }
                } // End of ORF loop

                // After processing all ORFs for this contig, write the buffer to the file stream.
                if (bb.length() > 0) {
                    bytesOutT += bb.length();
                    if (bsw != null) {
                        if (bsw.ordered) {
                            bsw.add(bb, r1.numericID);
                        } else {
                            bsw.addJob(bb);
                        }
                    }
                }
            } // End of contig loop

            // Finally, tell the system we are done with this batch of reads.
            cris.returnList(ln);
        }

        // Brandon - 7/23/25: Method to score an ORF with the neural network
        private float scoreOrf(Orf orf, Read contigRead) {
            fillFeatureVector(vec, orf, contigRead);
            net.applyInput(vec);
            return net.feedForward();
        }

        // Brandon - 7/23/25: Method to populate the feature vector for the network
        private void fillFeatureVector(float[] vec, Orf orf, Read contigRead) {
            Arrays.fill(vec, 0);
            int currentPos = 0;

            ContigStats contigStats = this.contigMetrics.get(contigRead.id);
            byte[] geneSeq = java.util.Arrays.copyOfRange(contigRead.bases, orf.start, orf.stop);
            if (orf.strand == 1) { AminoAcid.reverseComplementBasesInPlace(geneSeq); }

            // Feature 0: Scaled log length
            vec[currentPos++] = (float) (Math.log(orf.length()) / 10.0);
            // Feature 1: Start Score
            vec[currentPos++] = orf.startScore;
            // Feature 2: Avg Kmer Score
            vec[currentPos++] = orf.kmerScore / (orf.length() > 0 ? orf.length() : 1);
            // Feature 3: Stop Score
            vec[currentPos++] = orf.stopScore;
            // Feature 4: Gene GC
            vec[currentPos++] = CallGenesHelper.gcRatio(geneSeq, 0, geneSeq.length - 1);
            // Feature 5: Gene Entropy
            float geneEntropy = entropyTracker.averageEntropy(geneSeq, false);
            vec[currentPos++] = (float) Math.tanh(Math.log(geneEntropy / (1.0 + 1e-8 - geneEntropy)) / 4.0);
            // Feature 6: Contig GC
            vec[currentPos++] = (contigStats != null) ? (float) contigStats.gcRatio() : 0.0f;
            // Feature 7: Contig Entropy
            double contigEntropy = (contigStats != null) ? contigStats.entropy() : 0.0;
            vec[currentPos++] = (float) Math.tanh(Math.log(contigEntropy / (1.0 + 1e-8 - contigEntropy)) / 4.0);

            // One-hot encoded windows
            int length = orf.length();
            int mid = (orf.strand == 1) ? orf.start + (length / 2) - (length / 2 % 3) : orf.stop - (length / 2) + (length / 2 % 3);
            byte[] startWin = CallGenes.extractWindow(contigRead.bases, orf.start, 18, 18);
            byte[] midWin = CallGenes.extractWindow(contigRead.bases, mid, 6, 6);
            byte[] stopWin = CallGenes.extractWindow(contigRead.bases, orf.stop, 18, 18);

            if (orf.strand == 1) {
                AminoAcid.reverseComplementBasesInPlace(startWin);
                AminoAcid.reverseComplementBasesInPlace(midWin);
                AminoAcid.reverseComplementBasesInPlace(stopWin);
            }
            
            currentPos = fillOneHot(vec, currentPos, startWin);
            currentPos = fillOneHot(vec, currentPos, midWin);
            currentPos = fillOneHot(vec, currentPos, stopWin);
        }

        // Brandon - 7/23/25: Helper for filling one-hot vectors
        private int fillOneHot(float[] vec, int pos, byte[] bases) {
            for (byte b : bases) {
                char c = (char) b;
                switch (c) {
                    case 'A': case 'a': vec[pos] = 1; break;
                    case 'C': case 'c': vec[pos+1] = 1; break;
                    case 'G': case 'g': vec[pos+2] = 1; break;
                    case 'T': case 't': vec[pos+3] = 1; break;
                    default: break; // Ns are all zero
                }
                pos += 4;
            }
            return pos;
        }
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		ArrayList<Orf> processRead(final Read r){

			// Brandon:  Ensure all candidates are generated for vector output through this call (OLD CODE)
			/*if (vectorOut != null) {
                caller.setGenerateAllCandidates(true);
            }
			*/

			// Tell the GeneCaller to generate all candidates ONLY if the nofilter flag is set.
           	caller.setGenerateAllCandidates(this.nofilter);

		    // If nofilter=false, this list will be filtered by the caller's internal logic.
    		// If nofilter=true, this list will contain all potential candidates.
			ArrayList<Orf> list=caller.callGenes(r, pgm, true);

			// It's good practice to reset the flag after use 
			/*
            if (vectorOut != null) {
                caller.setGenerateAllCandidates(false);
            }
			*/
			
			if(geneHistT!=null && list!=null){
				for(Orf o : list){
					int bin=Tools.min(geneHistT.length-1, o.length()/geneHistDiv);
					geneHistT[bin]++;
				}
			}
			
			if(ros16S!=null){
				if(list!=null && !list.isEmpty()){
//					System.err.println("Looking for 16S.");
					ArrayList<Read> ssu=fetchType(r, list, r16S);
					if(ssu!=null && !ssu.isEmpty()){
//						System.err.println("Found "+ssu.size()+" 16S.");
						ros16S.add(ssu, r.numericID);
					}
				}
			}
			if(ros18S!=null){
				if(list!=null && !list.isEmpty()){
					ArrayList<Read> ssu=fetchType(r, list, r18S);
					if(ssu!=null && !ssu.isEmpty()){ros18S.add(ssu, r.numericID);}
				}
			}
			
			if(rosAmino!=null){
				if(mode==TRANSLATE){
					if(list!=null && !list.isEmpty()){
						ArrayList<Read> prots=translate(r, list);
						if(prots!=null){rosAmino.add(prots, r.numericID);}
					}
				}else if(mode==RETRANSLATE) {
					if(list!=null && !list.isEmpty()){
						ArrayList<Read> prots=translate(r, list);
						ArrayList<Read> ret=detranslate(prots);
						if(ret!=null){rosAmino.add(ret, r.numericID);}
					}
				}else if(mode==RECODE) {
					if(list!=null && !list.isEmpty()){
						Read recoded=recode(r, list);
						r.mate=null;
						ArrayList<Read> rec=new ArrayList<Read>(1);
						rec.add(recoded);
						if(rec!=null){rosAmino.add(rec, r.numericID);}
					}
				}else{
					assert(false) : mode;
				}
			}
			
			return list;
		}
		
		/** Number of reads processed by this thread */
		protected long readsInT=0;
		/** Number of bases processed by this thread */
		protected long basesInT=0;
		
		/** Number of genes called by this thread */
		protected long genesOutT=0;
		/** Number of bytes written by this thread */
		protected long bytesOutT=0;
		
		final long[] geneHistT;

		protected ConcurrentReadOutputStream rosAmino;
		protected ConcurrentReadOutputStream ros16S;
		protected ConcurrentReadOutputStream ros18S;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ByteStreamWriter bsw;
		/** Gene Model for annotation (not really needed) */
		private final GeneModel pgm;
		/** Gene Caller for annotation */
		final GeneCaller caller;
		/** Thread ID */
		final int tid;
		/** Entropy for this thread */
		final EntropyTracker entropyTracker;

		/** Flag to indicate if the header has been written for this thread's output */
		private boolean headerWritten = false;
		
		private int matchCount = 0; // Track number of Match rows for this thread
		
		public int getMatchCount() { return matchCount; }
	}
	
	public static ArrayList<Read> fetchType(final Read r, final ArrayList<Orf> list, int type){
		if(list==null || list.isEmpty()){return null;}
		ArrayList<Read> ret=new ArrayList<Read>(list.size());
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==type){
					Read sequence=fetch(orf, r.bases, r.id);
					ret.add(sequence);
				}
			}
			r.reverseComplement();
		}
		return (ret.size()>0 ? ret : null);
	}
	
	public static ArrayList<Read> translate(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return null;}
		ArrayList<Read> prots=new ArrayList<Read>(list.size());
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==CDS){
					Read aa=translate(orf, r.bases, r.id);
					prots.add(aa);
				}
			}
			r.reverseComplement();
		}
		return prots.isEmpty() ? null : prots;
	}
	
	public static Read recode(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return r;}
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==CDS){
					recode(orf, r.bases);
				}
			}
			r.reverseComplement();
		}
		return r;
	}
	
	public static ArrayList<Read> detranslate(final ArrayList<Read> prots){
		if(prots==null || prots.isEmpty()){return null;}
		ArrayList<Read> nucs=new ArrayList<Read>(prots.size());
		for(int strand=0; strand<2; strand++){
			for(Read prot : prots){
				Read nuc=detranslate(prot);
				nucs.add(nuc);
			}
		}
		return nucs;
	}
	
	public static Read translate(Orf orf, byte[] bases, String id){
//		assert(orf.length()%3==0) : orf.length(); //Happens sometimes on genes that go off the end, perhaps
		if(orf.strand==1){orf.flip();}
		byte[] acids=AminoAcid.toAAs(bases, orf.start, orf.stop);
		if(orf.strand==1){orf.flip();}
		Read r=new Read(acids, null, id+"\t"+(Shared.strandCodes[orf.strand]+"\t"+orf.start+"-"+orf.stop), 0, Read.AAMASK);
//		assert((r.length()+1)*3==orf.length());
		return r;
	}
	
	public static ArrayList<Read> fetch(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return null;}
		ArrayList<Read> genes=new ArrayList<Read>(list.size());
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==CDS){
					Read gene=fetch(orf, r.bases, r.id);
					genes.add(gene);
				}
			}
			r.reverseComplement();
		}
		return genes.isEmpty() ? null : genes;
	}
	
	public static Read fetch(Orf orf, Read source){
		assert(orf.start>=0 && orf.stop<source.length() && orf.length()>0) : 
			source.length()+"\n"+orf.length()+"\n"+orf;
		if(orf.strand==1){source.reverseComplement();}
		Read r=fetch(orf, source.bases, source.id);
		if(orf.strand==1){source.reverseComplement();}
		return r;
	}
	
	public static Read fetch(Orf orf, byte[] bases, String id){
		assert(orf.start>=0 && orf.stop<bases.length) : bases.length+"\n"+orf;
		if(orf.strand==1){orf.flip();}
		byte[] sub=Arrays.copyOfRange(bases, orf.start, orf.stop+1);
		if(orf.strand==1){orf.flip();}
		Read r=new Read(sub, null, id+"\t"+(Shared.strandCodes[orf.strand]+"\t"+orf.start+"-"+orf.stop), 0, 0);
		assert(r.length()==orf.length()) : r.length()+", "+orf.length();
		return r;
	}
	
	public static void recode(Orf orf, byte[] bases){
		if(orf.strand==1){orf.flip();}
		byte[] acids=AminoAcid.toAAs(bases, orf.start, orf.stop);
		for(int apos=0, bpos=orf.start; apos<acids.length; apos++){
			byte aa=acids[apos];
			int number=AminoAcid.acidToNumber[aa];
			String codon=(number>=0 ? AminoAcid.canonicalCodons[number] : "NNN");
			for(int i=0; i<3; i++, bpos++) {
				bases[bpos]=(byte)codon.charAt(i);
			}
		}
		if(orf.strand==1){orf.flip();}
	}
	
	public static Read detranslate(Read prot){
		ByteBuilder bb=new ByteBuilder(prot.length()*3);
		for(byte aa : prot.bases){
			int number=AminoAcid.acidToNumber[aa];
			String codon=(number>=0 ? AminoAcid.canonicalCodons[number] : "NNN");
			bb.append(codon);
		}
		Read r=new Read(bb.array, null, prot.id, prot.numericID, 0);
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static GeneCaller makeGeneCaller(GeneModel pgm){
		GeneCaller caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, 
				minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm);
		return caller;
	}
	
	private long maxReads=-1;
	private boolean merge;
	private boolean ecco;
	private int passes=1;
	
	private long readsIn=0;
	private long basesIn=0;
	private long genesOut=0;
	private long bytesOut=0;
	
	private static int minLen=80;//Actually a much higher value like 200 seems optimal compared to NCBI
	private static int maxOverlapSameStrand=80;
	private static int maxOverlapOppositeStrand=110;
	
	/* for kinnercds=6 */
//	private static float minStartScore=-0.10f;
//	private static float minStopScore=-0.5f;//Not useful; disabled
//	private static float minKmerScore=0.04f;//Does not seem useful.
//	private static float minOrfScore=40f; //Higher increases SNR dramatically but reduces TP rate
//	private static float minAvgScore=0.08f; //Not very effective

	/* for kinnercds=7 */
	private static float minStartScore=-0.10f;
	private static float minStopScore=-0.5f;//Not useful; disabled
	private static float minKmerScore=0.02f;//Does not seem useful.
	private static float minOrfScore=50f; //Higher increases SNR dramatically but reduces TP rate
	private static float minAvgScore=0.08f; //Not very effective
	
	long geneStopsMade=0;
	long geneStartsMade=0;
	long geneStartsRetained=0;
	long geneStopsRetained=0;
	long geneStartsOut=0;

	long tRNAOut=0;
	long r16SOut=0;
	long r23SOut=0;
	long r5SOut=0;
	long r18SOut=0;
	
	ScoreTracker stCds=new ScoreTracker(CDS);
	ScoreTracker stCds2=new ScoreTracker(CDS);
	ScoreTracker stCdsPass=new ScoreTracker(CDS);
	ScoreTracker sttRNA=new ScoreTracker(tRNA);
	ScoreTracker st16s=new ScoreTracker(r16S);
	ScoreTracker st23s=new ScoreTracker(r23S);
	ScoreTracker st5s=new ScoreTracker(r5S);
	ScoreTracker st18s=new ScoreTracker(r18S);
	
	ScoreTracker[] trackers=new ScoreTracker[] {stCdsPass, sttRNA, st16s, st23s, st5s, st18s};
	
	private int geneHistBins=1000;
	private int geneHistDiv=21;
	private final long[] geneHist;
	private boolean printZeroCountHist=false;
	
	/*--------------------------------------------------------------*/

	
	private HashSet<String> geneKeySet; // Set of gene keys loaded from true genes GFF file
	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> pgmList=new ArrayList<String>();
	private ArrayList<String> inGffList=new ArrayList<String>();
	private Map<String, ContigStats> contigMetrics; // Brandon: 6-24-25
	private HashSet<GeneQuad> trueGeneSet; // Brandon: 6-25-25
	private String outGff=null;
	private String outAmino=null;
	private String out16S=null;
	private String out18S=null;
	private String trueGenesFile = null;
	private String compareToGff=null;
	private String outStats="stderr";
	private String geneHistFile=null;
	private boolean json_out=false;
	// Brandon's CLI input calls
	private boolean seq = false;
	private boolean nofilter = false; 	// Brandon: 6/30/25

	//private boolean net = false; // Brandon: 7/23/25

    // Brandon - 7/23/25: Fields for neural network
    private String netFile = null;
    private final CellNet net0;
    private float cutoff = 0.5f;
    private boolean highpass = true;
	private FileFormat ffnet;

	private int totalGffRows = 0;
	private int totalGffGeneRows = 0;
	private int totalMatches = 0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffoutGff;
	private final FileFormat ffoutAmino;
	private final FileFormat ffout16S;
	private final FileFormat ffout18S;
	
	
	/** Determines how sequence is processed if it will be output */
	int mode=TRANSLATE;
	
	/** Translate nucleotides to amino acids */
	private static final int TRANSLATE=1;
	/** Translate nucleotides to amino acids,
	 * then translate back to nucleotides */
	private static final int RETRANSLATE=2;
	/** Re-encode coding regions of nucleotide
	 * sequences as a canonical codons */
	private static final int RECODE=3;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public boolean verbose=false;
	public boolean extendedStats=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean ordered=false; //this is OK sometimes, but sometimes hangs (e.g. on RefSeq mito), possibly if a sequence produces nothing.
	//To fix it, just ensure functions like translate always produce an array, even if it is empty.
	
	// Brandon's Functions
    /*
     * Convert a DNA sequence to a one-hot binary string (A=1000, C=0100, G=0010, T=0001, others=0000)
     */
	public static String toOneHotString(byte[] bases, int start, int end) {
		StringBuilder sb = new StringBuilder();
		for (int i = start; i <= end; i++) {
			char c = (char) bases[i];
			if (i > start) sb.append('\t');
			switch (c) {
				case 'A': case 'a':
					sb.append("1\t0\t0\t0");
					break;
				case 'C': case 'c':
					sb.append("0\t1\t0\t0");
					break;
				case 'G': case 'g':
					sb.append("0\t0\t1\t0");
					break;
				case 'T': case 't':
					sb.append("0\t0\t0\t1");
					break;
				default:
					sb.append("0\t0\t0\t0");
					break;
			}
		}
		return sb.toString();
	}

    /**
     * Extract a window of bases around a position, with padding as 'N' if out of bounds.
     * @param bases The base array
     * @param center The center position
     * @param before Number of bases before center (inclusive)
     * @param after Number of bases after center (inclusive)
     * @return byte[] window
     */
    public static byte[] extractWindow(byte[] bases, int center, int before, int after) {
        int len = before + after + 1;
        byte[] window = new byte[len];
        for (int i = 0; i < len; i++) {
            int idx = center - before + i;
            if (idx < 0 || idx >= bases.length) {
                window[i] = 'N';
            } else {
                window[i] = bases[idx];
            }
        }
        return window;
    }


	// Brandon
	// Convert a byte[] gene sequence to a String representation
	private static String basesToString(byte[] gene) {
	    StringBuilder sb = new StringBuilder();
	    for (byte base : gene) {
	        switch (base) {
	            case 'A': case 'a': sb.append('A'); break;
	            case 'C': case 'c': sb.append('C'); break;
	            case 'G': case 'g': sb.append('G'); break;
	            case 'T': case 't': sb.append('T'); break;
	            default: sb.append('N'); break;
	        }
	    }
	    return sb.toString();
	}



}
