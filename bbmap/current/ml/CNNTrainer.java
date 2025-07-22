package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import bin.Oracle;
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
import structures.ByteBuilder;
import structures.LongList;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Loads data for CNN training.
 * Based on ml.Trainer structure.
 * 
 * @author Brandon Imstepf 
 * Last Edited: @date 7-21-2025
 */
public class CNNTrainer {
    
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
        CNNTrainer x=new CNNTrainer(args);
        
        //Run the object
        x.process(t);
        
        //Close the print stream if it was redirected
        Shared.closeStream(x.outstream);
    }
    
    /**
     * Constructor.
     * @param args Command line arguments
     */
    public CNNTrainer(String[] args){
        
        {//Preparse block for help, config files, and outstream
            PreParser.printExecuting=false;
            Parser.printSetThreads=false;
            PreParser pp=new PreParser(args, getClass(), false);
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
            
            if(dataIn==null) {dataIn=parser.in1;}
        }
        
        validateParams();
        checkFileExistence(); //Ensure files can be read and written
        
        ffDataIn=FileFormat.testInput(dataIn, FileFormat.TXT, null, true, true);
        ffValidateIn=FileFormat.testInput(validateIn, FileFormat.TXT, null, true, true);
    }
    
    /*--------------------------------------------------------------*/
    /*----------------          Parsing             ----------------*/
    /*--------------------------------------------------------------*/
    
    Parser parse(String[] args) {
        
        Parser parser=new Parser();
        for(int i=0; i<args.length; i++){
            String arg=args[i];
            String[] split=arg.split("=");
            String a=split[0].toLowerCase();
            String b=split.length>1 ? split[1] : null;
            
            if(a.equals("data") || a.equals("train") || a.equals("traindata")){
                dataIn=b;
            }
            // CNN Architecture Parameters
            else if(a.equals("filters") || a.equals("numfilters")){
                // How many filters (feature detectors) in each conv layer
                // Example: filters=32,64,128
                String[] parts = b.split(",");
                filterCounts = new int[parts.length];
                for(int j=0; j<parts.length; j++){
                    filterCounts[j] = Integer.parseInt(parts[j]);
                }
            }else if(a.equals("filtersizes") || a.equals("kernels")){
                // Size of the sliding windows
                // Example: filtersizes=3,5,7
                String[] parts = b.split(",");
                filterSizes = new int[parts.length];
                for(int j=0; j<parts.length; j++){
                    filterSizes[j] = Integer.parseInt(parts[j]);
                }
            }else if(a.equals("poolsizes") || a.equals("pool")){
                // Pooling window sizes (usually 2)
                // Example: poolsizes=2,2,2
                String[] parts = b.split(",");
                poolSizes = new int[parts.length];
                for(int j=0; j<parts.length; j++){
                    poolSizes[j] = Integer.parseInt(parts[j]);
                }
            }else if(a.equals("dense") || a.equals("denselayers")){
                // Fully connected layer sizes
                // Example: dense=256,128
                String[] parts = b.split(",");
                denseLayers = new int[parts.length];
                for(int j=0; j<parts.length; j++){
                    denseLayers[j] = Integer.parseInt(parts[j]);
                }
            }
            // Training Parameters
            else if(a.equals("epochs")){
                epochs = Integer.parseInt(b);
            }else if(a.equals("batchsize") || a.equals("batch")){
                batchSize = Integer.parseInt(b);
            }else if(a.equals("learningrate") || a.equals("lr")){
                learningRate = Float.parseFloat(b);
            }else if(a.equals("dropout")){
                dropout = Float.parseFloat(b);
            }else if(a.equals("validate") || a.equals("validation") || a.equals("test")){
                validateIn=b;
            }else if(arg.equals("evaluate")){
                evaluate = true;
            }else if(a.equals("net") || a.equals("netin")){
                networkFile = b;
            }else if(a.equals("verbose")){
                verbose=Parse.parseBoolean(b);
            }else if(a.equals("shuffle")){
                shuffle=Parse.parseBoolean(b);
            }else if(a.equals("splitfraction") || a.equals("split")){
                splitFraction=Float.parseFloat(b);
            }else if(a.equals("maxlines")){
                maxLines=Parse.parseIntKMG(b);
            }else if(a.equals("vlines") || a.equals("vsamples")){
				maxLinesV=Parse.parseIntKMG(b);
				if(maxLinesV<0){maxLinesV=Integer.MAX_VALUE;}
			}else if(a.equals("exclusive")){
                exclusive=Parse.parseBoolean(b);
            }else if(a.equals("balance") || a.equals("balanced")){
				if(b==null || Tools.startsWithLetter(b)) {
					balance=Parse.parseBoolean(b) ? 1.0f : 0;
				}else{
					balance=Float.parseFloat(b);
				}
            }else if(parser.parse(arg, a, b)){
                //do nothing
            }else{
                outstream.println("Unknown parameter "+args[i]);
                assert(false) : "Unknown parameter "+args[i];
            }
        }

        
        
        return parser;
    }
    
    /*--------------------------------------------------------------*/
    /*----------------        Validation            ----------------*/
    /*--------------------------------------------------------------*/
    
    private void validateParams(){
        assert(splitFraction>=0 && splitFraction<=1) : "Split fraction must be between 0 and 1";
        if(maxLines<0){maxLines=Integer.MAX_VALUE;}
        if(maxLinesV<0){maxLinesV=maxLines;}
    }
    
    /*--------------------------------------------------------------*/
    /*----------------     File Existence           ----------------*/
    /*--------------------------------------------------------------*/
    
    private void checkFileExistence(){
        
        //Ensure input files can be read
        ArrayList<String> inFiles=new ArrayList<String>();
        if(dataIn!=null){inFiles.add(dataIn);}
        if(validateIn!=null){inFiles.add(validateIn);}
        
        if(inFiles.size()>0){
            if(!Tools.testInputFiles(false, true, inFiles.toArray(new String[0]))){
                throw new RuntimeException("\nCan't read some input files.\n");  
            }
        }
        
        //Ensure output files can be written
        ArrayList<String> outFiles=new ArrayList<String>();
        // Add output files here if needed
        
        if(!Tools.testOutputFiles(overwrite, append, false, outFiles.toArray(new String[0]))){
            throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files\n");
        }
    }
    
    /*--------------------------------------------------------------*/
    /*----------------         Outer Methods        ----------------*/
    /*--------------------------------------------------------------*/
    
    void process(Timer t){
    
    //Load training data
    SampleSet[] trainData = null;
    if(dataIn!=null){
        trainData = loadData(dataIn, true, maxLines, shuffle, splitFraction, maxLinesV);
        data = trainData[0];
        if(trainData.length>1 && validateSet==null){
            validateSet = trainData[1];
        }
    }
    
    //Load validation data if specified separately
    if(validateIn!=null){
        SampleSet[] valData = loadData(validateIn, false, maxLinesV, false, 0, 0);
        validateSet = valData[0];
    }
    
    //Validate that we have data
    if(data==null && validateSet==null){
        throw new RuntimeException("No data loaded!");
    }
    
    printStats(t);

    // Check if we're in evaluation mode
    if(evaluate) {
        // EVALUATION MODE
        outstream.println("\n" + "=".repeat(50));
        outstream.println("EVALUATION MODE");
        outstream.println("=".repeat(50));
        
        if(networkFile == null) {
            throw new RuntimeException("Evaluation mode requires a network file (net=filename)");
        }
        
        outstream.println("Loading network from: " + networkFile);
        CNNNetwork network = CNNNetwork.loadNetwork(networkFile);
        
        if(network == null) {
            throw new RuntimeException("Failed to load network from " + networkFile);
        }
        
        // Use validation set if available, otherwise use training data
        SampleSet evalData = (validateSet != null) ? validateSet : data;
        outstream.println("Evaluating on " + evalData.samples.length + " samples");
        
        // Run evaluation
        evaluateNetwork(network, evalData);
        
    } else {
        // TRAINING MODE (existing code)
        if(data != null) {
            outstream.println("\n" + "=".repeat(50));
            outstream.println("Initializing CNN Network...");
            
            // Create network with architecture parameters
            CNNNetwork network = new CNNNetwork(data.numInputs(), data.numOutputs());
            network.setArchitecture(filterCounts, filterSizes, poolSizes, denseLayers);
            network.setTrainingParams(epochs, batchSize, learningRate, dropout);
            network.initialize();
            
            outstream.println("\nStarting training...");
            network.train(data, validateSet);
            
            // Save network if output file specified
            if(networkOut != null) {
                outstream.println("\nSaving network to: " + networkOut);
                network.saveNetwork(networkOut);
            }
        }
    }
}
    
    /*--------------------------------------------------------------*/
    /*----------------         Inner Methods        ----------------*/
    /*--------------------------------------------------------------*/
    
    private SampleSet[] loadData(String path, boolean makeSubsets, int maxLines0, 
            boolean shuffleRaw, float splitFraction, int maxLines1){
        if(!quiet) {outstream.println("Loading "+path);}
        
        // Add debug output
        if(verbose) {
            outstream.println("Parameters: maxLines0="+maxLines0+", shuffle="+shuffleRaw+
                            ", split="+splitFraction+", maxLines1="+maxLines1+
                            ", exclusive="+exclusive+", balance="+balance);
        }
        
        SampleSet[] ssa = DataLoader.load(path, maxLines0, shuffleRaw, splitFraction, maxLines1, exclusive, balance);
        
        if(verbose) {
            outstream.println("DataLoader returned "+ssa.length+" sets");
        }
        
        for(SampleSet ss : ssa) {
            ss.makeSamples();
        }
        
        // Check if split actually worked
        if(ssa.length == 2 && splitFraction > 0) {
            int total = ssa[0].samples.length + ssa[1].samples.length;
            float actualSplit = (float)ssa[0].samples.length / total;
            
            if(verbose || Math.abs(actualSplit - splitFraction) > 0.1) {
                outstream.println("Split check: requested="+splitFraction+
                                ", actual="+String.format("%.3f", actualSplit)+
                                " ("+ssa[0].samples.length+"/"+total+")");
            }
        }
        
        // Original verbose output
        if(verbose) {
            for(int i=0; i<ssa.length; i++) {
                if(ssa[i] != null && ssa[i].samples != null) {
                    outstream.println("Set "+i+": "+ssa[i].samples.length+" samples, " +
                                    "positive="+ssa[i].numPositive+", negative="+ssa[i].numNegative);
                }
            }
        }
        
        if(makeSubsets && setsize>0 && ssa.length > 0) {
            assert(setsize>=100) : "Setsize should be at least 100";
            int subsets = Tools.max(1, ssa[0].samples.length/setsize);
            outstream.println("Data was organized into "+subsets+(subsets==1 ? " set." : " sets."));
            subsets = Tools.mid(1, subsets, ssa[0].samples.length);
            ssa[0].makeSubsets(subsets);
        }
        return ssa;
    }
        
    private void printStats(Timer t){
        outstream.println("\nData Loading Complete!");
        outstream.println("Time: \t"+t);
        
        if(data!=null){
            outstream.println("\nTraining Set:");
            outstream.println("Samples: \t"+data.samples.length);
            outstream.println("Inputs: \t"+data.numInputs());
            outstream.println("Outputs: \t"+data.numOutputs());
            outstream.println("Positive: \t"+data.numPositive);
            outstream.println("Negative: \t"+data.numNegative);
        }
        
        if(validateSet!=null){
            outstream.println("\nValidation Set:");
            outstream.println("Samples: \t"+validateSet.samples.length);
            outstream.println("Positive: \t"+validateSet.numPositive);
            outstream.println("Negative: \t"+validateSet.numNegative);
        }
        
        // Print first few samples if verbose
        if(verbose && data!=null && data.samples.length > 0){
            outstream.println("\nFirst 5 samples:");
            for(int i=0; i<Math.min(5, data.samples.length); i++){
                Sample s = data.samples[i];
                outstream.println("Sample "+i+": inputs="+s.in.length+
                                ", goal="+s.goal[0]+", positive="+s.positive);
            }
        }
    }

    private void evaluateNetwork(CNNNetwork network, SampleSet data) {
        Timer t = new Timer();
        
        int totalSamples = data.samples.length;
        int correct = 0;
        int tp = 0, fp = 0, tn = 0, fn = 0;
        float totalLoss = 0;
        
        // Arrays to store predictions and targets for ROC analysis
        float[] predictions = new float[totalSamples];
        float[] targets = new float[totalSamples];
        
        outstream.println("\nRunning evaluation...");
        
        // Process all samples
        for(int i = 0; i < totalSamples; i++) {
            Sample sample = data.samples[i];
            
            // Forward pass
            float[] output = network.forward(sample.in);
            float prediction = output[0];
            float target = sample.goal[0];
            
            predictions[i] = prediction;
            targets[i] = target;
            
            // Calculate loss
            totalLoss += BinaryCrossEntropyLoss.calculateLoss(prediction, target);
            
            // Binary classification metrics (using 0.5 threshold)
            boolean predicted = prediction > 0.5f;
            boolean actual = target > 0.5f;
            
            if(predicted == actual) correct++;
            
            if(actual && predicted) tp++;
            else if(!actual && !predicted) tn++;
            else if(!actual && predicted) fp++;
            else fn++;
        }
        
        // Calculate metrics
        float accuracy = (float)correct / totalSamples;
        float err = 1.0f - accuracy;
        float avgLoss = totalLoss / totalSamples;
        
        float precision = tp > 0 ? (float)tp / (tp + fp) : 0;
        float recall = tp > 0 ? (float)tp / (tp + fn) : 0;
        float f1 = (precision + recall) > 0 ? 
            2 * precision * recall / (precision + recall) : 0;
        
        float fpr = fp > 0 ? (float)fp / (fp + tn) : 0;
        float fnr = fn > 0 ? (float)fn / (fn + tp) : 0;
        float tpr = 1.0f - fnr;  // True positive rate (sensitivity)
        float tnr = 1.0f - fpr;  // True negative rate (specificity)
        
        // Calculate AUC (area under ROC curve) - simplified version
        float auc = calculateAUC(predictions, targets);
        
        // Output results
        outstream.println("\n==================================================");
        outstream.println("Evaluation Results:");
        outstream.println("==================================================");
        outstream.println(String.format("Samples:        %d", totalSamples));
        outstream.println(String.format("Positives:      %d (%.2f%%)", 
            tp + fn, 100.0f * (tp + fn) / totalSamples));
        outstream.println(String.format("Negatives:      %d (%.2f%%)", 
            tn + fp, 100.0f * (tn + fp) / totalSamples));
        outstream.println();
        
        outstream.println("Performance Metrics:");
        outstream.println(String.format("  Accuracy:     %.4f (%.2f%%)", accuracy, accuracy * 100));
        outstream.println(String.format("  Error Rate:   %.4f", err));
        outstream.println(String.format("  Avg Loss:     %.4f", avgLoss));
        outstream.println(String.format("  AUC:          %.4f", auc));
        outstream.println();
        
        outstream.println("Binary Classification (threshold=0.5):");
        outstream.println(String.format("  TP: %-6d  FP: %-6d", tp, fp));
        outstream.println(String.format("  TN: %-6d  FN: %-6d", tn, fn));
        outstream.println(String.format("  Precision:    %.4f", precision));
        outstream.println(String.format("  Recall (TPR): %.4f", recall));
        outstream.println(String.format("  F1 Score:     %.4f", f1));
        outstream.println(String.format("  FPR:          %.4f", fpr));
        outstream.println(String.format("  FNR:          %.4f", fnr));
        outstream.println(String.format("  TNR:          %.4f", tnr));
        
        // Threshold analysis
        outstream.println("\nThreshold Analysis:");
        float[] thresholds = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f};
        outstream.println("Threshold   FPR      FNR      F1       Accuracy");
        for(float thresh : thresholds) {
            MetricsAtThreshold m = calculateMetricsAtThreshold(predictions, targets, thresh);
            outstream.println(String.format("  %.1f     %.4f   %.4f   %.4f   %.4f", 
                thresh, m.fpr, m.fnr, m.f1, m.accuracy));
        }
        
        t.stop();
        outstream.println("\nEvaluation Time: " + t);
    }

    // Helper class for threshold metrics
    static class MetricsAtThreshold {
        float fpr, fnr, f1, accuracy;
    }

    // Calculate metrics at a specific threshold
    private MetricsAtThreshold calculateMetricsAtThreshold(
            float[] predictions, float[] targets, float threshold) {
        MetricsAtThreshold m = new MetricsAtThreshold();
        int tp = 0, fp = 0, tn = 0, fn = 0;
        
        for(int i = 0; i < predictions.length; i++) {
            boolean predicted = predictions[i] > threshold;
            boolean actual = targets[i] > 0.5f;
            
            if(actual && predicted) tp++;
            else if(!actual && !predicted) tn++;
            else if(!actual && predicted) fp++;
            else fn++;
        }
        
        m.fpr = fp > 0 ? (float)fp / (fp + tn) : 0;
        m.fnr = fn > 0 ? (float)fn / (fn + tp) : 0;
        m.accuracy = (float)(tp + tn) / predictions.length;
        
        float precision = tp > 0 ? (float)tp / (tp + fp) : 0;
        float recall = tp > 0 ? (float)tp / (tp + fn) : 0;
        m.f1 = (precision + recall) > 0 ? 
            2 * precision * recall / (precision + recall) : 0;
        
        return m;
    }

    // Simple AUC calculation
    private float calculateAUC(float[] predictions, float[] targets) {
        // Count positives and negatives
        int positives = 0, negatives = 0;
        for(float t : targets) {
            if(t > 0.5f) positives++;
            else negatives++;
        }
        
        if(positives == 0 || negatives == 0) return 0.5f;
        
        // Simple Mann-Whitney U statistic approach
        long sum = 0;
        for(int i = 0; i < predictions.length; i++) {
            if(targets[i] > 0.5f) { // positive sample
                for(int j = 0; j < predictions.length; j++) {
                    if(targets[j] < 0.5f) { // negative sample
                        if(predictions[i] > predictions[j]) sum += 1;
                        else if(predictions[i] == predictions[j]) sum += 0.5;
                    }
                }
            }
        }
        
        return (float)sum / (positives * negatives);
    }
    
    /*--------------------------------------------------------------*/
    /*----------------            Fields            ----------------*/
    /*--------------------------------------------------------------*/
    
    /** Primary input file path */
    private String dataIn=null;
    
    /** Validation file path */
    private String validateIn=null;
    
    /** Primary output stream */
    private PrintStream outstream=System.err;
    
    /** Maximum lines to load from training data */
    private int maxLines=-1;
    
    /** Maximum lines to load from validation data */
    private int maxLinesV=-1;
    
    /** Shuffle the data */
    private boolean shuffle=true;
    
    /** Train/validation split fraction */
    private float splitFraction=0.1f;
    
    /** Set size for organizing data into subsets */
    private int setsize=1000;
    
    /** Exclusive mode */
    private boolean exclusive=true;
    
    /** Balance positive/negative samples */
    private float balance=0.0f;
    
    /** Verbose output */
    private boolean verbose=false;
    
    /** Don't print typical output */
    private boolean quiet=false;
    
    /** Overwrite existing output files */
    private boolean overwrite=false;
    
    /** Append to existing output files */
    private boolean append=false;

    /*  Evaluate the network on a validation dataset, rather than train */
    boolean evaluate = false;

    /* Existing network file to load */
    String networkFile = null;

    private String networkOut = null; // Output file for saving the trained network

    /** CNN Architecture Parameters */
    private int[] filterCounts = {32, 64};  // Number of filters per conv layer
    private int[] filterSizes = {5, 3};     // Size of convolution windows
    private int[] poolSizes = {2, 2};       // Pooling window sizes
    private int[] denseLayers = {128};      // Fully connected layer sizes

    /** Training Parameters */
    private int epochs = 50;                 // How many times to go through all data
    private int batchSize = 32;             // Samples per mini-batch
    private float learningRate = 0.001f;    // How fast to learn
    private float dropout = 0.5f;           // Randomly drop neurons to prevent overfitting

    
    /*--------------------------------------------------------------*/
    /*----------------        Data Structures       ----------------*/
    /*--------------------------------------------------------------*/
    
    /** Training data */
    private SampleSet data;
    
    /** Validation data */
    private SampleSet validateSet;
    
    /** File format for data input */
    private FileFormat ffDataIn;
    
    /** File format for validation input */
    private FileFormat ffValidateIn;
}