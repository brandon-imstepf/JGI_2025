/* Brandon Imstepf
* @date 6/24/25
* */

package prok;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Reader;
import java.io.InputStream;
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
import java.util.Set;

import prok.Orf;
import stream.Read;
import dna.AminoAcid;
import tracker.EntropyTracker;


public class CallGenesHelper{

    public record ContigStats(double gcRatio, double entropy) {}

    /*
    * @date 6-24-25
    * Reads a FASTA file into a map of contig IDs to sequences.
    * This version correctly handles both plain-text and gzipped FASTA files.
    * @param filePath Path to the FASTA file.
    * @return A map where keys are contig IDs and values are sequences.
    */
    public static Map<String, String> readFastaFile(String filePath) throws IOException {
        Map<String, String> sequences = new HashMap<>();
        BufferedReader reader;

        // Check if the file is gzipped and choose the correct way to open it
        if (filePath.endsWith(".gz")) {
            InputStream fileStream = new FileInputStream(filePath);
            InputStream gzipStream = new GZIPInputStream(fileStream); // Decompress the stream
            Reader decoder = new InputStreamReader(gzipStream);
            reader = new BufferedReader(decoder);
        } else {
            reader = new BufferedReader(new FileReader(filePath));
        }

        // This try-with-resources block ensures the reader (whether plain or gzipped) is closed.
        try (reader) {
            String line;
            String currentContigId = null;
            StringBuilder currentSequence = new StringBuilder();

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (currentContigId != null) {
                        sequences.put(currentContigId, currentSequence.toString());
                    }
                    currentContigId = line.substring(1).trim();
                    currentSequence.setLength(0);
                } else {
                    currentSequence.append(line.trim());
                }
            }
            // Store the last contig in the file
            if (currentContigId != null) {
                sequences.put(currentContigId, currentSequence.toString());
            }
        }
        return sequences;
    }


    /**
     * @date 6-24-25
     * Pre-calculates statistics for all contigs.
     * @param contigSequences The map of contig IDs to sequences.
     * @return A map of contig IDs to their calculated stats.
     */
    public static Map<String, ContigStats> calculateContigMetrics(Map<String, String> contigSequences) {
        Map<String, ContigStats> contigMetrics = new HashMap<>();
        for (Map.Entry<String, String> entry : contigSequences.entrySet()) {
            String contigId = entry.getKey();
            String sequence = entry.getValue();

            double gcRatio = calculateGCRatio(sequence); 

            // Example: k=5, window=50, nucleotide mode (amino=false)
            EntropyTracker et = new EntropyTracker(5, 50, false);
            double entropy = et.averageEntropy(sequence.getBytes(), false);

            contigMetrics.put(contigId, new ContigStats(gcRatio, entropy));
        }
        System.err.println("Pre-calculated statistics for " + contigMetrics.size() + " contigs.");
        return contigMetrics;
    }

    /**
     * @date 6-24-25
     * Calculates the entropy of a DNA sequence.
     * @param dna The DNA sequence.
     * @return The entropy value.
     */
    private static double calculateGCRatio(String dna) {
        if (dna == null || dna.isEmpty()) { return 0.0; }
        long gcCount = 0;
        long atcgCount = 0;
        for (int i = 0; i < dna.length(); i++) {
            char c = Character.toUpperCase(dna.charAt(i));
            if (c == 'G' || c == 'C') {
                gcCount++;
                atcgCount++;
            } else if (c == 'A' || c == 'T') {
                atcgCount++;
            }
        }
        return (atcgCount > 0) ? (double) gcCount / atcgCount : 0.0;
    }

    /** Brandon
     * @date 6-18-25
     * Calculates GC-content for a DNA sequence string.
     * @param bases The DNA sequence as a byte array.
     * @param from The starting index (inclusive).
     * @param to The ending index (inclusive).
     * @return The GC-ratio from 0.0 to 1.0.
     */
    public static float gcRatio(byte[] bases, int from, int to) {
        int gc = 0, len = 0;
        for (int i = from; i <= to && i < bases.length; i++) {
            byte b = bases[i];
            if (b == 'G' || b == 'g' || b == 'C' || b == 'c') { gc++; }
            if (b == 'A' || b == 'a' || b == 'T' || b == 't' || b == 'G' || b == 'g' || b == 'C' || b == 'c') { len++; }
        }
        return len > 0 ? (float)gc / len : 0f;
    }


    /**
     * Generates a feature vector string for a given ORF.
     * This centralizes all feature calculation and formatting logic.
     *
     * @param orf The ORF (gene candidate) to process.
     * @param contigRead The full Read object for the contig the ORF belongs to.
     * @param contigStats The pre-calculated statistics for the entire contig.
     * @param entropyTracker An instance of EntropyTracker to calculate gene entropy.
     * @param geneKeySet A set of "true gene" coordinates for ground truth matching. Can be null.
     * @param seqMode If true, output the raw DNA sequence instead of one-hot encoding.
     * @param mlMode If true, do not include the contig ID prefix.
     * @return A tab-separated string of features for the ORF.
     */
    public static String generateFeatureVector(
            Orf orf,
            Read contigRead,
            ContigStats contigStats,
            EntropyTracker entropyTracker,
            Set<String> geneKeySet,
            boolean seqMode,
            boolean mlMode) {

        // --- 1. Get Gene Sequence and Basic Features ---
        int start = orf.start;
        int stop = orf.stop;
        String strand = (orf.strand == 1 ? "-" : "+");
        byte[] geneSeq = java.util.Arrays.copyOfRange(contigRead.bases, start, stop); // Use stop, as copyOfRange is exclusive

        // The full contig read is always on the forward strand, so we only flip the gene sequence itself if needed.
        if (orf.strand == 1) { 
            AminoAcid.reverseComplementBasesInPlace(geneSeq);
        }
        
        // --- 2. Calculate Gene-Specific Metrics ---
        float geneGc = gcRatio(geneSeq, 0, geneSeq.length - 1);
        float geneEntropy = entropyTracker.averageEntropy(geneSeq, false);
        
        // --- 3. Get Contig-Specific Metrics ---
        double contigGc = (contigStats != null) ? contigStats.gcRatio() : 0.0;
        double contigEntropy = (contigStats != null) ? contigStats.entropy() : 0.0;
        
        // --- 4. Ground Truth Matching ---
        String match = "0";
        if (geneKeySet != null) {
            String key = (start + 1) + "\t" + (stop) + "\t" + strand; // Adjusted to match your GFF loading logic
            if (geneKeySet.contains(key)) {
                match = "1";
            } 
            // Note: Fuzzy goes here for later
        }
        
        // --- 5. Build the Feature String ---
        StringBuilder sb = new StringBuilder();
        if (!mlMode) {
            sb.append(contigRead.id).append('\t');
        }
        
        // Debug: Add the start and stop coordinates to compare:
        //sb.append(start + 1).append('\t'); 
        //sb.append(stop+1).append('\t'); 
        // Append all the calculated features in order.
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", Math.log(orf.length()) / 10.0)).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", orf.startScore)).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", orf.kmerScore / (orf.length() > 0 ? orf.length() : 1))).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", orf.stopScore)).append('\t');
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", geneGc)).append('\t');
        //sb.append(String.format(java.util.Locale.ROOT, "%.6f", geneEntropy)).append('\t');
        double geneEntropyLogitScaled = Math.tanh(Math.log(geneEntropy / (1.0 + 1e-8 - geneEntropy)) / 4.0);
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", geneEntropyLogitScaled)).append('\t');
        // Add the new contig-level stats!
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", contigGc)).append('\t');
        //sb.append(String.format(java.util.Locale.ROOT, "%.6f", contigEntropy));
        double contigEntropyLogitScaled = Math.tanh(Math.log(contigEntropy / (1.0 + 1e-8 - contigEntropy)) / 4.0);
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", contigEntropyLogitScaled)).append('\t');

        // --- 6. Handle Sequence/One-Hot Encoding ---
        if (seqMode) {
            sb.append('\t').append(new String(geneSeq));
        } else {
            // This logic is from your original code.
            int length = orf.length();
            int mid = (orf.strand == 1) ? start + (length / 2) - (length / 2 % 3) : stop - (length / 2) + (length / 2 % 3);
            
            // Note: The window extraction needs the original contig bases and coordinates
            byte[] startWin = CallGenes.extractWindow(contigRead.bases, start, 18, 18);
            byte[] midWin = CallGenes.extractWindow(contigRead.bases, mid, 6, 6);
            byte[] stopWin = CallGenes.extractWindow(contigRead.bases, stop, 18, 18);

            if (orf.strand == 1) {
                AminoAcid.reverseComplementBasesInPlace(startWin);
                AminoAcid.reverseComplementBasesInPlace(midWin);
                AminoAcid.reverseComplementBasesInPlace(stopWin);
            }

            sb.append('\t').append(CallGenes.toOneHotString(startWin, 0, startWin.length - 1));
            sb.append('\t').append(CallGenes.toOneHotString(midWin, 0, midWin.length - 1));
            sb.append('\t').append(CallGenes.toOneHotString(stopWin, 0, stopWin.length - 1));
        }
        
        // --- 7. Append Match Status ---
        if (geneKeySet != null) {
            sb.append('\t').append(match);
        }
        
        return sb.toString();
    }



}