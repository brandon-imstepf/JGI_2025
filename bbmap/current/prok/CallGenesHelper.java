package prok;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.Collections;

import dna.AminoAcid;
import ml.CellNet;
import ml.CellNetParser;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import fileIO.ByteStreamWriter;
import tracker.EntropyTracker;
import java.io.File;
import java.util.List;

import shared.Shared;
import gff.GffLine;

public class CallGenesHelper {

    // RECORDS
    public record ContigStats(double gcRatio, double entropy) {}
    public record GeneQuad(String contig, int start, int stop, byte strand) {}
    private record TrueGeneData(HashSet<GeneQuad> set, int totalRows, int geneRows) {}
    private record ScoreData(String orfId, float originalScore, float modifiedScore) {}
    public record GffStats(int truePositives, int falsePositives, int falseNegatives, int refCount) {}

    // MEMBER FIELDS
    private String trueGenesFile = null;
    private String netFile = null;
    private String compareToGff = null;
    private boolean nofilter = false;
    private boolean seqMode = false;
    private CellNet net0;
    private Map<String, ContigStats> contigMetrics;
    private HashSet<GeneQuad> trueGeneSet;
    private int totalGffRows = 0;
    private int totalGffGeneRows = 0;
    private int totalMatches = 0;
    private CellNet threadNet;
    private float[] threadVec;
    private EntropyTracker threadEntropyTracker;
    private int threadMatchCount = 0;
    private float cutoff = 0.5f;
    private boolean highpass = true;
    
    // Neural Network Fields
    private static CellNet primaryNet = null;
    private static float nnCutoff = 0.5f;
    private static float nnStrength = 0.5f;
    private static boolean nnDebug = false;
    private static float nnMinScore = 0.0f;
    private static boolean nnAllOrfs = false;
    private static float prefilterCutoff = 0.1f;
    
    // Debug and assertion counters
    private static int netLoadAssertions = 0;
    private static int threadCloneAssertions = 0;
    private static int scoreModificationCalls = 0;
    
    // Thread-local neural network fields
    private CellNet threadNeuralNet = null;
    private float[] threadFeatureVector = null;

    // Logging Fields
    private boolean enableLogging = true;
    private String logDirectoryPath = null;
    private final ArrayList<Float> finalScores = new ArrayList<>();
    private int truePositives = 0;
    private int falsePositives = 0;
    private int threadTruePositives = 0;
    private int threadFalsePositives = 0;
    private ByteStreamWriter scoreCsvStream = null;
    private Set<GeneQuad> calledCdsQuads = new HashSet<>();
    private Set<GeneQuad> threadCalledCdsQuads = new HashSet<>();

    // PUBLIC METHODS

    public boolean isLoggingEnabled() {
        return enableLogging;
    }

    public String generateScoresCsvPath(String outGff) {
        if (!enableLogging || outGff == null) {
            return null;
        }
        File logDir = getLogDirectory(outGff);
        if (logDir == null) { return null; }
        return new File(logDir, "scores.csv").getPath();
    }

    public boolean parse(String arg, String a, String b) {
        if (a.equals("log")) {
            enableLogging = shared.Parse.parseBoolean(b);
            return true;
        }
        if (a.equals("truegenes")) { trueGenesFile = b; return true; }
    // Accept both forms: 'nofilter' (flag) and 'nofilter=true' (key=value)
    if (a.equalsIgnoreCase("nofilter") || arg.equalsIgnoreCase("nofilter")) { nofilter = true; return true; }
        if (a.equals("net")) { netFile = b; return true; }
        if (arg.equalsIgnoreCase("seq")) { seqMode = true; return true; }
        if (a.equals("compareto")) { compareToGff = b; return true; }
        
        // Neural Network Parameters
        if (a.equals("net") || a.equals("neuralnet") || a.equals("nn")) {
            netFile = b;
            System.err.println("DEBUG: Neural network file set to: " + netFile);
            return true;
        }
        if (a.equals("cutoff") || a.equals("nncutoff")) {
            nnCutoff = Float.parseFloat(b);
            System.err.println("DEBUG: Neural network cutoff set to: " + nnCutoff);
            assert(nnCutoff > 0 && nnCutoff < 1) : "Neural network cutoff must be between 0 and 1, got: " + nnCutoff;
            return true;
        }
        if (a.equals("strength") || a.equals("nnstrength")) {
            nnStrength = Float.parseFloat(b);
            System.err.println("DEBUG: Neural network strength set to: " + nnStrength);
            assert(nnStrength >= 0 && nnStrength <= 1) : "Neural network strength must be between 0 and 1, got: " + nnStrength;
            return true;
        }
        if (a.equals("nndebug")) {
            nnDebug = shared.Parse.parseBoolean(b);
            System.err.println("DEBUG: Neural network debug mode: " + nnDebug);
            return true;
        }
        if (a.equals("nnminscore")) {
            nnMinScore = Float.parseFloat(b);
            return true;
        }
        if (a.equals("nnallorfs")) {
            nnAllOrfs = shared.Parse.parseBoolean(b);
            return true;
        }
        if (a.equals("prefilter_cutoff")) {
            prefilterCutoff = Float.parseFloat(b);
            return true;
        }
        
        return false;
    }

    public void initialize(ArrayList<String> fnaList, PrintStream outstream, String outGff, String compareToGff) {
        this.compareToGff = compareToGff;
        if (trueGenesFile != null) {
            try {
                TrueGeneData data = loadTrueGeneSet(trueGenesFile);
                this.trueGeneSet = data.set();
                this.totalGffRows = data.totalRows();
                this.totalGffGeneRows = data.geneRows();
            } catch (IOException e) { outstream.println("ERROR: Failed to load true genes file: " + trueGenesFile); e.printStackTrace(); }
        }
        if (netFile != null) {
            net0 = CellNetParser.load(netFile);
            assert(net0 != null) : "Net file is null or incorrectly loaded: " + netFile;
        }
        
        // Load neural network if specified
        if (netFile != null) {
            loadNeuralNetwork(outstream);
        }
        
        if (fnaList != null && !fnaList.isEmpty()) {
            try {
                String fastaFile = fnaList.get(0);
                Map<String, String> contigSequences = readFastaFile(fastaFile);
                this.contigMetrics = calculateContigMetrics(contigSequences);
            } catch (IOException e) { outstream.println("ERROR: Failed to read FASTA for contig metrics."); e.printStackTrace(); this.contigMetrics = new HashMap<>(); }
        }

        if (enableLogging && outGff != null) {
            File logDir = getLogDirectory(outGff);
            if (!logDir.exists()) {
                logDir.mkdirs();
            }
            this.logDirectoryPath = logDir.getAbsolutePath();
        }
    }
    
    public CallGenesHelper getThreadLocalCopy() {
        CallGenesHelper copy = new CallGenesHelper();
        copy.net0 = this.net0;
        copy.contigMetrics = this.contigMetrics;
        assert(this.contigMetrics == null || copy.contigMetrics != null) : "contigMetrics was lost during thread copy.";
        copy.trueGeneSet = this.trueGeneSet;
        copy.nofilter = this.nofilter;
        copy.seqMode = this.seqMode;
        copy.cutoff = this.cutoff;
        copy.highpass = this.highpass;
        copy.trueGenesFile = this.trueGenesFile;
        copy.enableLogging = this.enableLogging;
        copy.compareToGff = this.compareToGff;
        copy.calledCdsQuads = this.calledCdsQuads;
        
        // Clone neural network for thread safety
        if (primaryNet != null) {
            // ASSERTION: Ensure primary network exists before cloning
            assert(primaryNet != null) : "THREAD_CLONE_ASSERTION: Primary neural network is null during thread copy";
            
            // Clone the neural network for thread safety
            copy.threadNeuralNet = primaryNet.copy(false);
            
            // ASSERTION: Ensure cloning was successful
            assert(copy.threadNeuralNet != null) : "THREAD_CLONE_ASSERTION: Failed to clone neural network for thread";
            assert(copy.threadNeuralNet != primaryNet) : "THREAD_CLONE_ASSERTION: Cloned network is same object as primary";
            assert(copy.threadNeuralNet.numInputs() == primaryNet.numInputs()) : "THREAD_CLONE_ASSERTION: Input dimensions mismatch";
            assert(copy.threadNeuralNet.numOutputs() == primaryNet.numOutputs()) : "THREAD_CLONE_ASSERTION: Output dimensions mismatch";
            
            // Initialize feature vector based solely on the network's input size.
            int netInputs = copy.threadNeuralNet.numInputs();
            copy.threadFeatureVector = new float[netInputs];
            if (nnDebug) {
                System.err.println("DEBUG: Allocated feature vector of size " + netInputs);
                System.err.println("DEBUG: Thread successfully cloned neural network");
            }           
            // Increment assertion counter
            threadCloneAssertions++;
        }
        
        return copy;
    }

    public void initializeThreadObjects() {
        if (net0 != null) {
            this.threadNet = net0.copy(false);
            this.threadVec = new float[net0.numInputs()];
        }
        this.threadEntropyTracker = new EntropyTracker(5, 50, false);
    }
    

   public ArrayList<Orf> processContig(Read contig, GeneCaller caller, GeneModel pgm,
                                     ConcurrentReadOutputStream rosAmino, ConcurrentReadOutputStream ros16S, ConcurrentReadOutputStream ros18S) {
     
        final String originalId = contig.id;
        if (contig.id != null && contig.id.contains(" ")) {
            contig.id = contig.id.split("\\s+")[0];
        }
        assert contig.id != null && !contig.id.contains(" ") : "Failed to trim contig ID: " + originalId;

        // Set up neural network integration in GeneCaller
        caller.setHelper(this);
        caller.setCurrentContigRead(contig);

        ArrayList<Orf> finalOrfs;
        // Only generate all candidates if 'nofilter' is true
        if (nofilter) {
            // RADAR: Assert that we are intentionally entering a special processing mode.
            assert(nofilter) : "Logic error: should not be in 'all candidates' mode without nofilter flag.";

            caller.setGenerateAllCandidates(true);
            ArrayList<Orf> candidates = caller.callGenes(contig, pgm, true);
            caller.setGenerateAllCandidates(false); // Reset

            // When nofilter is true, return all candidates without filtering
            finalOrfs = candidates;
        } else { // Default Mode (with or without neural network score modification)
            finalOrfs = caller.callGenes(contig, pgm, true);
        }

        if (finalOrfs != null && !finalOrfs.isEmpty()) {
            handleSequenceOutputs(contig, finalOrfs, rosAmino, ros16S, ros18S);

            if (enableLogging) {
                for (Orf orf : finalOrfs) {
                    if (orf.type == ProkObject.CDS) {
                        this.finalScores.add(orf.orfScore); // Add to thread-local list
                        if (trueGeneSet != null) {
                            GeneQuad orfQuad = new GeneQuad(contig.id, orf.start + 1, orf.stop + 1, (byte)orf.strand);
                            if (trueGeneSet.contains(orfQuad)) {
                                threadTruePositives++;
                            } else {
                                threadFalsePositives++;
                            }
                        }
                        threadCalledCdsQuads.add(new GeneQuad(contig.id, orf.start + 1, orf.stop + 1, (byte)orf.strand));
                    }
                }
            }
        }
        
        return finalOrfs;
    }

    public synchronized void accumulateStatsFromThread(CallGenesHelper threadHelper) {
        this.totalMatches += threadHelper.threadMatchCount;
        if(this.enableLogging) {
            this.truePositives += threadHelper.threadTruePositives;
            this.falsePositives += threadHelper.threadFalsePositives;
            this.finalScores.addAll(threadHelper.finalScores);
            this.calledCdsQuads.addAll(threadHelper.threadCalledCdsQuads);
        }
    }
    
    public void printFinalStats(PrintStream outstream) {
        if (enableLogging) {
            writeLogFile();
            if (scoreCsvStream != null && trueGeneSet != null) {
                appendFalseNegativesToScores(scoreCsvStream);
            }
        }
        if (trueGenesFile != null) {
            outstream.println();
            outstream.println("True Genes File: " + trueGenesFile);
            outstream.println("Total Rows detected in .gff = " + totalGffRows);
            outstream.println(".gff Gene Rows = " + totalGffGeneRows);
            outstream.println("Total Matches = " + totalMatches);
        }
    }

    public void formatGffOutput(ByteBuilder bb, ArrayList<Orf> orfList, Read contig) {
        assert(bb != null) : "ByteBuilder cannot be null.";
        assert(orfList != null) : "Orf list cannot be null.";
        assert(contig != null && contig.id != null) : "Contig and its ID cannot be null.";

        final boolean trainingMode = (this.trueGenesFile != null && this.net0 == null);
        final String contigId = contig.id.split("\\s+")[0];

        for (Orf orf : orfList) {
            assert(orf != null) : "ORF cannot be null.";
            // When in training mode, we are only interested in labeling and
            // outputting the protein-coding genes (CDS).
            if (trainingMode && orf.type != ProkObject.CDS) {
                continue; // Skip tRNAs, rRNAs, etc.
            }

            // Manually build the GFF line using the correct API calls.
            bb.append(contigId).tab().append("BBTools").tab().append(ProkObject.typeStrings2[orf.type]).tab();
            bb.append(orf.start + 1).tab().append(orf.stop + 1).tab();
            bb.append(orf.orfScore, 2).tab();
            bb.append(orf.strand < 0 ? '.' : Shared.strandCodes2[orf.strand]).tab();
            bb.append(orf.frame).tab();
            
            // Manually build the rich attribute string.
            bb.append("ID=").append(contigId).append('_').append(orf.start).append(';');
            bb.append(ProkObject.typeStrings[orf.type]);
            if (orf.type == ProkObject.CDS) { bb.append(",fr").append(orf.frame); }
            bb.append(",startScr:").append(orf.startScore, 3);
            bb.append(",stopScr:").append(orf.stopScore, 3);
            bb.append(",innerScr:").append(orf.averageKmerScore(), 3);
            bb.append(",len:").append(orf.length());
            if (orf.type == ProkObject.CDS) {
                bb.append(',').append("start:").append(AminoAcid.codonToString(orf.startCodon));
                bb.append(',').append("stop:").append(AminoAcid.codonToString(orf.stopCodon));
            }

            // Append the VECTOR attribute. This 'if' is now redundant but safe to keep.
            if (trainingMode && orf.type == ProkObject.CDS) {
                ContigStats stats = contigMetrics.get(contigId);
                assert(stats != null) : "ContigStats missing for key: " + contigId;
                String vectorString = generateFeatureVector(orf, contig, stats, threadEntropyTracker, trueGeneSet, seqMode);
                if (vectorString.endsWith("\t1")) { this.threadMatchCount++; }
                bb.append(";VECTOR=");
                for (int i = 0; i < vectorString.length(); i++) {
                    char c = vectorString.charAt(i);
                    bb.append(c == '\t' ? ',' : c);
                }
            }
            bb.nl();
        }
    }

    private void writeLogFile() {
        if (logDirectoryPath == null) {
            return;
        }

        // Don't create a log file if no genes were scored.
        if (finalScores.isEmpty() && trueGeneSet == null && compareToGff == null) {
            return;
        }

        try (PrintStream ps = new PrintStream(new File(logDirectoryPath, "log.txt"))) {
            
            // FP/FN rates
            if (trueGeneSet != null && !trueGeneSet.isEmpty()) {
                int totalCalledCds = truePositives + falsePositives;
                float fpRate = (totalCalledCds > 0) ? (float)falsePositives / totalCalledCds : 0;
                
                int totalTrueGenes = totalGffGeneRows;
                int falseNegatives = totalTrueGenes - truePositives;
                float fnRate = (totalTrueGenes > 0) ? (float)falseNegatives / totalTrueGenes : 0;

                ps.println("False positive rate: " + String.format(java.util.Locale.ROOT, "%.6f", fpRate));
                ps.println("False negative rate: " + String.format(java.util.Locale.ROOT, "%.6f", fnRate));
                ps.println("True positives: " + truePositives);
                ps.println("False positives: " + falsePositives);
                ps.println("False negatives: " + falseNegatives);
            } else if (compareToGff != null) {
                try {
                    String outGffPath = new File(logDirectoryPath).getParent() + "/" + new File(logDirectoryPath).getName() + ".gff";
                    GffStats stats = gradeGff(outGffPath, compareToGff);
                    
                    float fpRate = (stats.truePositives() + stats.falsePositives() > 0) ? (float)stats.falsePositives() / (stats.truePositives() + stats.falsePositives()) : 0;
                    float fnRate = (stats.refCount() > 0) ? (float)stats.falseNegatives() / stats.refCount() : 0;

                    ps.println("False positive rate: " + String.format(java.util.Locale.ROOT, "%.6f", fpRate));
                    ps.println("False negative rate: " + String.format(java.util.Locale.ROOT, "%.6f", fnRate));
                    ps.println("True positives: " + stats.truePositives());
                    ps.println("False positives: " + stats.falsePositives());
                    ps.println("False negatives: " + stats.falseNegatives());
                } catch (IOException e) {
                    ps.println("Could not generate grading stats: " + e.getMessage());
                }
            }
            
            // Cutoff
            ps.println("Cutoff: " + nnCutoff);
            
            // Scores
            if(!finalScores.isEmpty()){
                Collections.sort(finalScores);
                float minScore = finalScores.get(0);
                float maxScore = finalScores.get(finalScores.size() - 1);
                
                double sum = 0;
                for (float score : finalScores) {
                    sum += score;
                }
                double avgScore = sum / finalScores.size();
                
                float medianScore;
                int size = finalScores.size();
                if (size % 2 == 0) {
                    medianScore = (finalScores.get(size / 2 - 1) + finalScores.get(size / 2)) / 2.0f;
                } else {
                    medianScore = finalScores.get(size / 2);
                }

                ps.println("Highest score: " + maxScore);
                ps.println("Lowest score: " + minScore);
                ps.println("Average score: " + avgScore);
                ps.println("Median score: " + medianScore);
            }

        } catch (IOException e) {
            // It's a log file, so if it fails, just print to stderr.
            System.err.println("Failed to write log file.");
            e.printStackTrace();
        }
    }

    private GffStats gradeGff(String queryGffPath, String refGffPath) throws IOException {
        Map<GeneQuad, GffLine> refMap = loadGffToMap(refGffPath);
        Map<GeneQuad, GffLine> queryMap = loadGffToMap(queryGffPath);

        int truePositives = 0;
        int falsePositives = 0;

        for (GeneQuad queryQuad : queryMap.keySet()) {
            if (refMap.containsKey(queryQuad)) {
                truePositives++;
            } else {
                falsePositives++;
            }
        }

        int falseNegatives = 0;
        for (GeneQuad refQuad : refMap.keySet()) {
            if (!queryMap.containsKey(refQuad)) {
                falseNegatives++;
            }
        }

        return new GffStats(truePositives, falsePositives, falseNegatives, refMap.size());
    }

    public void setScoreWriter(ByteStreamWriter bsw) {
        this.scoreCsvStream = bsw;
    }

    /**
     * Converts an integer ORF type from ProkObject to its GFF string representation.
     * This is a self-contained alternative to relying on a potentially private
     * API in the ProkObject class.
     * @param type The integer type of the ORF (e.g., ProkObject.CDS).
     * @return The string representation (e.g., "CDS").
     */
    private static String orfTypeToString(int type) {
        // RADAR: Asserting the precondition that the type is valid.
        assert(type >= 0) : "Invalid ORF type: " + type;

        switch (type) {
            case ProkObject.CDS: return "CDS";
            case ProkObject.tRNA: return "tRNA";
            case ProkObject.r16S: return "rRNA"; // Common practice to label all rRNAs as 'rRNA'
            case ProkObject.r18S: return "rRNA";
            case ProkObject.r23S: return "rRNA";
            case ProkObject.r5S: return "rRNA";
            default:
                // This assertion acts as a firewall. If a new, unhandled type is ever
                // added, the program will crash here instead of producing a corrupt GFF.
                assert(false) : "Unhandled ORF type: " + type;
                return "unknown";
        }
    }
    
    private String generateFeatureVector(Orf orf, Read contigRead, ContigStats contigStats, 
                                     EntropyTracker entropyTracker, Set<GeneQuad> trueGeneSet, boolean seqMode) {
                                        
        assert(contigStats != null) : 
                "ContigStats lookup failed for contig ID: '" + contigRead.id + "'. " +
                "Attempted lookup with short key: '" + contigRead.id.split("\\s+")[0] + "'. " +
                "This key was not found in the metrics map. Available keys: " + contigMetrics.keySet();
    
        int start = orf.start;
        int stop = orf.stop;
        String match = "0";

        if (trueGeneSet != null) {
            byte orfStrand = (byte) orf.strand;
            String contigIdShort = contigRead.id.split("\\s+")[0];
            GeneQuad orfQuad = new GeneQuad(contigIdShort, start + 1, stop + 1, orfStrand);
            if (trueGeneSet.contains(orfQuad)) {
                match = "1";
            }
        }

        byte[] geneSeq = java.util.Arrays.copyOfRange(contigRead.bases, start, stop + 1);
        if (orf.strand == 1) { AminoAcid.reverseComplementBasesInPlace(geneSeq); }

        // These calculations will now be correct
        float geneGc = gcRatio(geneSeq, 0, geneSeq.length - 1);
        float geneEntropy = entropyTracker.averageEntropy(geneSeq, false);
        double contigGc = (contigStats != null) ? contigStats.gcRatio() : 0.0;
        double contigEntropy = (contigStats != null) ? contigStats.entropy() : 0.0;
        
        StringBuilder sb = new StringBuilder();
        sb.append(contigRead.id).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", Math.log(orf.length()) / 10.0)).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", orf.startScore)).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", orf.kmerScore / (orf.length() > 0 ? orf.length() : 1))).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", orf.stopScore)).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", geneGc)).append('\t');
        double geneEntropyLogitScaled = Math.tanh(Math.log(geneEntropy / (1.0 + 1e-8 - geneEntropy)) / 4.0);
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", geneEntropyLogitScaled)).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", contigGc)).append('\t');
        double contigEntropyLogitScaled = Math.tanh(Math.log(contigEntropy / (1.0 + 1e-8 - contigEntropy)) / 4.0);
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", contigEntropyLogitScaled));

        if (seqMode) {
            sb.append('\t').append(new String(geneSeq));
        } else {
            int length = orf.length();
            int offset = (length / 2) % 3;
            int mid = (orf.strand == 0)
                ? start + (length / 2 - offset)
                : stop - (length / 2 - offset);            
            byte[] startWin = extractWindow(contigRead.bases, start, 18, 18);
            byte[] midWin = extractWindow(contigRead.bases, mid, 6, 6); // The missing window
            byte[] stopWin = extractWindow(contigRead.bases, stop, 18, 18);

            if (orf.strand == 1) {
                AminoAcid.reverseComplementBasesInPlace(startWin);
                AminoAcid.reverseComplementBasesInPlace(midWin);
                AminoAcid.reverseComplementBasesInPlace(stopWin);
            }
            sb.append('\t').append(toOneHotString(startWin));
            sb.append('\t').append(toOneHotString(midWin)); // Add it to the output
            sb.append('\t').append(toOneHotString(stopWin));
        }

        if (trueGeneSet != null) { sb.append('\t').append(match); }
        return sb.toString();
    }
    
    private void handleSequenceOutputs(Read r, ArrayList<Orf> list, ConcurrentReadOutputStream rosAmino, ConcurrentReadOutputStream ros16S, ConcurrentReadOutputStream ros18S) {
        // This logic is safe to keep static calls as it uses original CallGenes methods
        if (ros16S != null) { ArrayList<Read> ssu = CallGenes.fetchType(r, list, ProkObject.r16S); if(ssu != null && !ssu.isEmpty()){ ros16S.add(ssu, r.numericID); } }
        if (ros18S != null) { ArrayList<Read> ssu = CallGenes.fetchType(r, list, ProkObject.r18S); if(ssu != null && !ssu.isEmpty()){ ros18S.add(ssu, r.numericID); } }
        if (rosAmino != null) { ArrayList<Read> prots = CallGenes.translate(r, list); if(prots != null){ rosAmino.add(prots, r.numericID); } }
    }
    
    private float scoreOrf(Orf orf, Read contigRead) {
        // Use neural network if available
        if (hasNeuralNetwork()) {
            modifyOrfScoreWithNeuralNetwork(orf, contigRead);
            return orf.orfScore;
        }
        // Fallback to original score
        return orf.orfScore;
    }
    
    // --- STATIC HELPER METHODS ---
    
private static Map<String, String> readFastaFile(String filePath) throws IOException {
    Map<String, String> sequences = new HashMap<>();
    try (BufferedReader reader = openBufferedReader(filePath)) {
        String line;
        String currentContigId = null;
        StringBuilder currentSequence = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if (currentContigId != null) { 
                    sequences.put(currentContigId, currentSequence.toString()); 
                }
                // Ensure we only use the ID before the first whitespace.
                currentContigId = line.substring(1).split("\\s+")[0];
                currentSequence.setLength(0);
            } else {
                currentSequence.append(line.trim());
            }
        }
        if (currentContigId != null) { 
            sequences.put(currentContigId, currentSequence.toString()); 
        }
    }
    return sequences;
}

    private static Map<String, ContigStats> calculateContigMetrics(Map<String, String> contigSequences) {
        Map<String, ContigStats> contigMetrics = new HashMap<>();
        EntropyTracker et = new EntropyTracker(5, 50, false);
        for (Map.Entry<String, String> entry : contigSequences.entrySet()) {
            String sequence = entry.getValue();
            double gcRatio = gcRatio(sequence.getBytes(), 0, sequence.length() -1);
            double entropy = et.averageEntropy(sequence.getBytes(), false);
            contigMetrics.put(entry.getKey(), new ContigStats(gcRatio, entropy));
        }
        return contigMetrics;
    }

    private static TrueGeneData loadTrueGeneSet(String gffFilePath) throws IOException {
        int totalRows = 0, geneRows = 0;
        HashSet<GeneQuad> trueGeneSet = new HashSet<>();
        try (BufferedReader br = openBufferedReader(gffFilePath)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.trim().isEmpty()) continue;
                totalRows++;
                String[] cols = line.split("\t");
                if (cols.length > 7 && cols[2].equalsIgnoreCase("CDS")) {
                    geneRows++;
                    try {
                        String contig = cols[0];
                        int start = Integer.parseInt(cols[3]);
                        int stop = Integer.parseInt(cols[4]);
                        byte strand = (byte) (cols[6].equals("+") ? 0 : 1);
                        trueGeneSet.add(new GeneQuad(contig, start, stop, strand));
                    } catch (NumberFormatException e) { /* Malformed line, ignore */ }
                }
            }
        }
        return new TrueGeneData(trueGeneSet, totalRows, geneRows);
    }
    
    private static BufferedReader openBufferedReader(String filename) throws IOException {
        if (filename.endsWith(".gz")) {
            return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename))));
        } else {
            return new BufferedReader(new FileReader(filename));
        }
    }

    private static float gcRatio(byte[] bases, int from, int to) {
        int gc = 0, len = 0;
        for (int i = from; i <= to && i < bases.length; i++) {
            byte b = bases[i];
            if (b == 'G' || b == 'g' || b == 'C' || b == 'c') { gc++; }
            if (b == 'A' || b == 'a' || b == 'T' || b == 't' || b == 'G' || b == 'g' || b == 'C' || b == 'c') { len++; }
        }
        return len > 0 ? (float)gc / len : 0f;
    }

    private static byte[] extractWindow(byte[] bases, int center, int before, int after) {
        int len = before + after + 1;
        byte[] window = new byte[len];
        for (int i = 0; i < len; i++) {
            int idx = center - before + i;
            if (idx < 0 || idx >= bases.length) { window[i] = 'N'; }
            else { window[i] = bases[idx]; }
        }
        return window;
    }

    private static String toOneHotString(byte[] bases) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < bases.length; i++) {
            char c = (char) bases[i];
            if (i > 0) sb.append('\t');
            switch (c) {
                case 'A': case 'a': sb.append("1\t0\t0\t0"); break;
                case 'C': case 'c': sb.append("0\t1\t0\t0"); break;
                case 'G': case 'g': sb.append("0\t0\t1\t0"); break;
                case 'T': case 't': sb.append("0\t0\t0\t1"); break;
                default: sb.append("0\t0\t0\t0"); break;
            }
        }
        return sb.toString();
    }
    
    // --- NEURAL NETWORK METHODS ---
    
    /**
     * Load the neural network from file with comprehensive assertions
     */
    private void loadNeuralNetwork(PrintStream outstream) {
        // Check if neural network file is specified
        if (netFile == null || netFile.trim().isEmpty()) {
            if (nnDebug) {
                outstream.println("DEBUG: No neural network file specified, skipping neural network loading");
            }
            return;
        }
        
        try {
            // Load the neural network
            primaryNet = CellNetParser.load(netFile);
            
            // Check if loading was successful
            if (primaryNet == null) {
                outstream.println("WARNING: Failed to load neural network from file: " + netFile);
                return;
            }
            
            if (primaryNet.numInputs() <= 0 || primaryNet.numOutputs() <= 0) {
                outstream.println("WARNING: Neural network has invalid dimensions - inputs: " +
                    primaryNet.numInputs() + ", outputs: " + primaryNet.numOutputs());
                primaryNet = null;
                return;
            }
            
            // Debug output
            if (nnDebug) {
                outstream.println("DEBUG: Successfully loaded neural network from: " + netFile);
                outstream.println("DEBUG: Network inputs: " + primaryNet.numInputs());
                outstream.println("DEBUG: Network outputs: " + primaryNet.numOutputs());
                outstream.println("DEBUG: Neural network cutoff: " + nnCutoff);
                outstream.println("DEBUG: Neural network strength: " + nnStrength);
            }
            
            // Increment assertion counter
            netLoadAssertions++;
            
        } catch (Exception e) {
            outstream.println("WARNING: Failed to load neural network from: " + netFile + " - " + e.getMessage());
            if (nnDebug) {
                e.printStackTrace();
            }
            primaryNet = null;
        }
    }
    
    /**
     * Modify ORF score using neural network prediction
     * Algorithm: mult = NN_output / cutoff; final_mult = ((mult - 1) × strength) + 1; modified_score = original_score × final_mult
     */
    public void modifyOrfScoreWithNeuralNetwork(Orf orf, Read contigRead) {
        // ASSERTION: Ensure we have a neural network available
        assert(threadNeuralNet != null) : "SCORE_MOD_ASSERTION: Thread neural network is null";
        assert(threadFeatureVector != null) : "SCORE_MOD_ASSERTION: Thread feature vector is null";
        assert(orf != null) : "SCORE_MOD_ASSERTION: ORF is null";
        assert(contigRead != null) : "SCORE_MOD_ASSERTION: Contig read is null";
        
        final float originalScore = orf.orfScore;

        // Phase 2: Prefilter based on score/length ratio
        if (!nnAllOrfs && orf.length() > 0 && (originalScore / orf.length()) < prefilterCutoff) {
            orf.orfScore = -1; // Explicitly mark the ORF as invalid
            return;
        }
        
        // Increment call counter first
        scoreModificationCalls++;
        
        if (!nnAllOrfs && originalScore < nnMinScore) {
            return;
        }
        
        // Debug output every 10000 calls to reduce spam but still track progress
        if (scoreModificationCalls % 10000 == 0) {
            System.err.println("DEBUG: Neural network called " + scoreModificationCalls + " times. Current ORF: " + orf.start + "-" + orf.stop + ", score: " + originalScore);
        }
        
        try {
            // OPTIMIZED: Generate feature vector directly into float array to avoid string operations
            generateFeatureVectorDirect(orf, contigRead, threadFeatureVector);
            
            // Run neural network inference
            threadNeuralNet.applyInput(threadFeatureVector);
            threadNeuralNet.feedForward();
            float[] output = threadNeuralNet.getOutput();
            
            // ASSERTION: Ensure output is valid
            assert(output != null && output.length > 0) : "SCORE_MOD_ASSERTION: Neural network output is null or empty";
            assert(!Float.isNaN(output[0]) && !Float.isInfinite(output[0])) : "SCORE_MOD_ASSERTION: Neural network output is NaN or infinite: " + output[0];
            
            // Apply score modification algorithm
            float mult = output[0] / nnCutoff;
            float finalMult = ((mult - 1) * nnStrength) + 1;
            float modifiedScore = originalScore * finalMult;
            
            // ASSERTION: Ensure modified score is reasonable
            assert(!Float.isNaN(modifiedScore) && !Float.isInfinite(modifiedScore)) : "SCORE_MOD_ASSERTION: Modified score is NaN or infinite: " + modifiedScore;
            
            // Update the ORF's score directly
            orf.orfScore = modifiedScore;
            
            // Debug output for first few calls only
            if (nnDebug && scoreModificationCalls <= 5) {
                System.err.println("DEBUG: NN output: " + output[0] + ", mult: " + mult + ", finalMult: " + finalMult);
                System.err.println("DEBUG: Original score: " + originalScore + " -> Modified score: " + modifiedScore);
            }
            
            // Log score data for CSV
            if (scoreCsvStream != null) {
                String orfId = contigRead.id + "_" + orf.start;
                
                String status = "Unknown";
                if (trueGeneSet != null) {
                    GeneQuad orfQuad = new GeneQuad(contigRead.id, orf.start + 1, orf.stop + 1, (byte)orf.strand);
                    if (trueGeneSet.contains(orfQuad)) {
                        status = "TP";
                    } else {
                        status = "FP";
                    }
                }

                ByteBuilder bb = new ByteBuilder();
                bb.append(orfId).append(',');
                bb.append(String.format(java.util.Locale.ROOT, "%.4f", originalScore)).append(',');
                bb.append(String.format(java.util.Locale.ROOT, "%.4f", modifiedScore)).append(',');
                bb.append(orf.length()).append(',');
                bb.append(status).nl();
                scoreCsvStream.add(bb, contigRead.numericID);
               }
            
        } catch (Exception e) {
            if (nnDebug) {
                System.err.println("ERROR: Neural network score modification failed: " + e.getMessage());
                e.printStackTrace();
            }
            
            // ASSERTION: Neural network operations should not fail
            assert(false) : "SCORE_MOD_ASSERTION: Exception during neural network score modification: " + e.getMessage();
        }
    }
    
    /**
     * Optimized feature vector generation that directly populates float array
     * Avoids expensive string operations and parsing
     */
    private void generateFeatureVectorDirect(Orf orf, Read contigRead, float[] featureVector) {
        // Get cached contig stats
        String contigIdShort = contigRead.id.split("\\s+")[0];
        ContigStats stats = contigMetrics.get(contigIdShort);
        assert(stats != null) : "ContigStats missing for contig: " + contigIdShort;
        
        int start = orf.start;
        int stop = orf.stop;
        
        // Extract gene sequence (this is still needed for GC and entropy calculations)
        byte[] geneSeq = java.util.Arrays.copyOfRange(contigRead.bases, start, stop + 1);
        if (orf.strand == 1) {
            AminoAcid.reverseComplementBasesInPlace(geneSeq);
        }
        
        // Calculate features directly into array (skip contig ID which was index 0 in string version)
        int idx = 0;
        
        // Feature 1: log(length) / 10.0
        if (idx < featureVector.length) featureVector[idx++] = (float)(Math.log(orf.length()) / 10.0);
        
        // Feature 2: start score
        if (idx < featureVector.length) featureVector[idx++] = orf.startScore;
        
        // Feature 3: average kmer score
        if (idx < featureVector.length) featureVector[idx++] = orf.kmerScore / (orf.length() > 0 ? orf.length() : 1);
        
        // Feature 4: stop score
        if (idx < featureVector.length) featureVector[idx++] = orf.stopScore;
        
        // Feature 5: gene GC ratio
        if (idx < featureVector.length) featureVector[idx++] = gcRatio(geneSeq, 0, geneSeq.length - 1);
        
        // Feature 6: gene entropy (logit scaled)
        if (idx < featureVector.length) {
            float geneEntropy = threadEntropyTracker.averageEntropy(geneSeq, false);
            double geneEntropyLogitScaled = Math.tanh(Math.log(geneEntropy / (1.0 + 1e-8 - geneEntropy)) / 4.0);
            featureVector[idx++] = (float)geneEntropyLogitScaled;
        }
        
        // Feature 7: contig GC ratio
        if (idx < featureVector.length) featureVector[idx++] = (float)stats.gcRatio();
        
        // Feature 8: contig entropy (logit scaled)
        if (idx < featureVector.length) {
            double contigEntropy = stats.entropy();
            double contigEntropyLogitScaled = Math.tanh(Math.log(contigEntropy / (1.0 + 1e-8 - contigEntropy)) / 4.0);
            featureVector[idx++] = (float)contigEntropyLogitScaled;
        }
        
        // Features 9+: One-hot encoded sequence windows (if not in seq mode)
        // This is the most expensive part - we need to optimize this
        int length = orf.length();
        int mid = (orf.strand == 0) ? start + (length / 2) - (length / 2 % 3) : stop - (length / 2) + (length / 2 % 3);
        
        // Extract windows
        byte[] startWin = extractWindow(contigRead.bases, start, 18, 18);
        byte[] midWin = extractWindow(contigRead.bases, mid, 6, 6);
        byte[] stopWin = extractWindow(contigRead.bases, stop, 18, 18);
        
        if (orf.strand == 1) {
            AminoAcid.reverseComplementBasesInPlace(startWin);
            AminoAcid.reverseComplementBasesInPlace(midWin);
            AminoAcid.reverseComplementBasesInPlace(stopWin);
        }
        
        // Convert to one-hot directly into feature vector with bounds checking
        idx = encodeOneHotDirect(startWin, featureVector, idx);
        idx = encodeOneHotDirect(midWin, featureVector, idx);
        idx = encodeOneHotDirect(stopWin, featureVector, idx);
        
        // Debug: Check if we exceeded bounds
        if (idx > featureVector.length) {
            System.err.println("WARNING: Feature vector overflow. Expected: " + featureVector.length + ", got: " + idx);
        }
    }
    
    /**
     * Optimized one-hot encoding directly into float array
     */
    private int encodeOneHotDirect(byte[] bases, float[] featureVector, int startIdx) {
        int idx = startIdx;
        for (byte base : bases) {
            // Ensure there is enough space for 4 floats per base
            if (idx + 4 > featureVector.length) {
                if (nnDebug) {
                    System.err.println("WARNING: Feature vector overflow while encoding base at index " + idx);
                }
                break;
            }
            char c = (char) base;
            switch (c) {
                case 'A': case 'a':
                    featureVector[idx++] = 1.0f; featureVector[idx++] = 0.0f;
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 0.0f;
                    break;
                case 'C': case 'c':
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 1.0f;
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 0.0f;
                    break;
                case 'G': case 'g':
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 0.0f;
                    featureVector[idx++] = 1.0f; featureVector[idx++] = 0.0f;
                    break;
                case 'T': case 't':
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 0.0f;
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 1.0f;
                    break;
                default:
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 0.0f;
                    featureVector[idx++] = 0.0f; featureVector[idx++] = 0.0f;
                    break;
            }
        }
        return idx;
    }

    /**
     * Check if neural network is available for this thread
     */
    public boolean hasNeuralNetwork() {
        return threadNeuralNet != null;
    }
    
    /**
     * Get neural network status information
     */
    public String getNeuralNetworkStatus() {
        if (primaryNet == null) {
            return "No neural network loaded";
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("Neural Network Status:\n");
        sb.append("  Primary network loaded: ").append(primaryNet != null).append("\n");
        sb.append("  Thread network available: ").append(threadNeuralNet != null).append("\n");
        sb.append("  Network file: ").append(netFile != null ? netFile : "none").append("\n");
        sb.append("  Cutoff: ").append(nnCutoff).append("\n");
        sb.append("  Strength: ").append(nnStrength).append("\n");
        sb.append("  Debug mode: ").append(nnDebug).append("\n");
        sb.append("  Load assertions: ").append(netLoadAssertions).append("\n");
        sb.append("  Clone assertions: ").append(threadCloneAssertions).append("\n");
        sb.append("  Score modification calls: ").append(scoreModificationCalls);
        
        return sb.toString();
    }
    
    /**
     * Reset neural network debug counters
     */
    public static void resetNeuralNetworkCounters() {
        netLoadAssertions = 0;
        threadCloneAssertions = 0;
        scoreModificationCalls = 0;
    }

    private static Map<GeneQuad, GffLine> loadGffToMap(String filePath) throws IOException {
        Map<GeneQuad, GffLine> map = new HashMap<>();
        try (BufferedReader reader = openBufferedReader(filePath)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                
                GffLine gffLine = new GffLine(line.getBytes());

                if ("CDS".equalsIgnoreCase(gffLine.type)) {
                    try {
                        GeneQuad quad = new GeneQuad(gffLine.seqid, gffLine.start, gffLine.stop, (byte)gffLine.strand);
                        map.put(quad, gffLine);
                    } catch (Exception e) {
                        System.err.println("Warning: Skipping malformed GFF line: " + line);
                    }
                }
            }
        }
        return map;
    }

    // --- PRIVATE HELPER METHODS ---

    private File getLogDirectory(String outGff) {
        if (outGff == null) { return null; }
        File gffFile = new File(outGff);
        String name = gffFile.getName();
        if (name.endsWith(".gz")) { name = name.substring(0, name.length() - 3); }
        if (name.endsWith(".gff")) { name = name.substring(0, name.length() - 4); }
        if (name.endsWith(".tmp")) { name = name.substring(0, name.length() - 4); }
        File parentDir = gffFile.getParentFile();
        return (parentDir == null) ? new File(name) : new File(parentDir, name);
    }

    private void appendFalseNegativesToScores(ByteStreamWriter writer) {
        if (trueGeneSet == null || calledCdsQuads == null) {
            return;
        }

        // Find genes in true set that were not in the called set
        Set<GeneQuad> falseNegatives = new HashSet<>(trueGeneSet);
        falseNegatives.removeAll(calledCdsQuads);

        for (GeneQuad fn : falseNegatives) {
            String orfId = fn.contig() + "_" + (fn.start() - 1);
            int length = fn.stop() - fn.start() + 1;
            ByteBuilder bb = new ByteBuilder();
            bb.append(orfId).append(',');
            bb.append("0.0000").append(','); // Original Score
            bb.append("0.0000").append(','); // Modified Score
            bb.append(length).append(','); // Length
            bb.append("FN").nl(); // Status
            writer.add(bb, 0); // Use a dummy numeric ID
        }
    }
}