package prok;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import dna.AminoAcid;
import prok.Orf;
import stream.Read;
import tracker.EntropyTracker;

public class CallGenesHelper {

    public record ContigStats(double gcRatio, double entropy) {}
    public record GeneQuad(String contig, int start, int stop, byte strand) {}
    public record TrueGeneData(HashSet<GeneQuad> set, int totalRows, int geneRows) {}

    public static Map<String, String> readFastaFile(String filePath) throws IOException {
        Map<String, String> sequences = new HashMap<>();
        BufferedReader reader;

        if (filePath.endsWith(".gz")) {
            InputStream fileStream = new FileInputStream(filePath);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream);
            reader = new BufferedReader(decoder);
        } else {
            reader = new BufferedReader(new FileReader(filePath));
        }

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
            if (currentContigId != null) {
                sequences.put(currentContigId, currentSequence.toString());
            }
        }
        return sequences;
    }

    public static Map<String, ContigStats> calculateContigMetrics(Map<String, String> contigSequences) {
        Map<String, ContigStats> contigMetrics = new HashMap<>();
        EntropyTracker et = new EntropyTracker(5, 50, false);
        for (Map.Entry<String, String> entry : contigSequences.entrySet()) {
            String contigId = entry.getKey();
            String sequence = entry.getValue();
            double gcRatio = calculateGCRatio(sequence);
            double entropy = et.averageEntropy(sequence.getBytes(), false);
            contigMetrics.put(contigId, new ContigStats(gcRatio, entropy));
        }
        System.err.println("Pre-calculated statistics for " + contigMetrics.size() + " contigs.");
        return contigMetrics;
    }

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

    public static float gcRatio(byte[] bases, int from, int to) {
        int gc = 0, len = 0;
        for (int i = from; i <= to && i < bases.length; i++) {
            byte b = bases[i];
            if (b == 'G' || b == 'g' || b == 'C' || b == 'c') { gc++; }
            if (b == 'A' || b == 'a' || b == 'T' || b == 't' || b == 'G' || b == 'g' || b == 'C' || b == 'c') { len++; }
        }
        return len > 0 ? (float)gc / len : 0f;
    }

    public static TrueGeneData loadTrueGeneSet(String gffFilePath) throws IOException {
        int totalRows = 0;
        int geneRows = 0;
        HashSet<GeneQuad> trueGeneSet = new HashSet<>();

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gffFilePath))))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.trim().isEmpty()) continue;
                totalRows++;
                String[] cols = line.split("\t");
                if (cols.length > 7 && cols[2].equalsIgnoreCase("CDS")) {
                    geneRows++;
                    try {
                        // *** THIS IS THE FIX ***
                        // Apply the same logic here as in generateFeatureVector to ensure keys match.
                        String contig = cols[0].split("\\s+")[0];
                        
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

    public static String generateFeatureVector(
            Orf orf,
            Read contigRead,
            ContigStats contigStats,
            EntropyTracker entropyTracker,
            Set<GeneQuad> trueGeneSet,
            boolean seqMode,
            boolean mlMode) {

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

        double geneEntropyLogitScaled = Math.tanh(Math.log(geneEntropy / (1.0 + 1e-8 - geneEntropy)) / 4.0);
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", geneEntropyLogitScaled)).append('\t');

        sb.append(String.format(java.util.Locale.ROOT, "%.6f", contigGc)).append('\t');
        double contigEntropyLogitScaled = Math.tanh(Math.log(contigEntropy / (1.0 + 1e-8 - contigEntropy)) / 4.0);
        sb.append(String.format(java.util.Locale.ROOT, "%.6f", contigEntropyLogitScaled));

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
