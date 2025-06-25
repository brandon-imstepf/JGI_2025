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
    // Hold coordinates of a gene in a contig.
    public record GeneQuad(String contig, int start, int stop, byte strand) {}
    // A container to return multiple values from our GFF loading method.
    public record TrueGeneData(HashSet<GeneQuad> set, int totalRows, int geneRows) {}

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
     * Loads a GFF file of "true" genes into a HashSet of GeneQuad objects for fast lookups.
     * This entire process now lives inside the helper class.
     * @param gffFilePath Path to the gzipped GFF file.
     * @return A TrueGeneData object containing the HashSet and row counts.
     */
    public static TrueGeneData loadTrueGeneSet(String gffFilePath) throws IOException {
        int totalRows = 0;
        int geneRows = 0;
        HashSet<GeneQuad> trueGeneSet = new HashSet<>();

        // Using the smart gzipped reader we built earlier
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gffFilePath))))) {
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
                    } catch (NumberFormatException e) {
                        System.err.println("Warning: Could not parse coordinates in GFF line: " + line);
                    }
                }
            }
        }
        return new TrueGeneData(trueGeneSet, totalRows, geneRows);
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
    // In CallGenesHelper.java

    public static String generateFeatureVector(
            Orf orf,
            Read contigRead,
            ContigStats contigStats,
            EntropyTracker entropyTracker,
            Set<GeneQuad> trueGeneSet, // This parameter is correctly named here
            boolean seqMode,
            boolean mlMode) {

        int start = orf.start;
        int stop = orf.stop;

        String match = "0";
        if (trueGeneSet != null) { // We are checking the parameter 'trueGeneSet'
            byte orfStrand = (byte) orf.strand;
            // The GFF start is 1-based, the Orf start is 0-based.
            GeneQuad orfQuad = new GeneQuad(contigRead.id, start + 1, stop, orfStrand);
            if (trueGeneSet.contains(orfQuad)) { // We are using the parameter 'trueGeneSet'
                match = "1";
            }
        }

        // --- The rest of the function remains the same ---
        byte[] geneSeq = java.util.Arrays.copyOfRange(contigRead.bases, start, stop);
        if (orf.strand == 1) { AminoAcid.reverseComplementBasesInPlace(geneSeq); }
        float geneGc = gcRatio(geneSeq, 0, geneSeq.length - 1);
        float geneEntropy = entropyTracker.averageEntropy(geneSeq, false);
        double contigGc = (contigStats != null) ? contigStats.gcRatio() : 0.0;
        double contigEntropy = (contigStats != null) ? contigStats.entropy() : 0.0;
        StringBuilder sb = new StringBuilder();
        if (!mlMode) { sb.append(contigRead.id).append('\t'); }
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
        
        if (seqMode) {
            sb.append('\t').append(new String(geneSeq));
        } else {
            int length = orf.length();
            int mid = (orf.strand == 1) ? start + (length / 2) - (length / 2 % 3) : stop - (length / 2) + (length / 2 % 3);
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
        if (trueGeneSet != null) { sb.append('\t').append(match); }
        return sb.toString();
    }



}