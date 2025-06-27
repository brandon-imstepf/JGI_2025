package prok;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader; // FIX: Added missing import
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import gff.GffLine;
import prok.CallGenesHelper.GeneQuad; // Reusing our nifty GeneQuad record

/**
 * A command-line tool to perform set operations on GFF files based on feature coordinates.
 * @author The GFF Engineering Team (Brandon Imstepf, et al.)
 * @date June 25, 2025
 */
public class GffSetOperator {

    public static void main(String[] args) {
        // --- Argument Parsing ---
        String fileA = null;
        String fileB = null;
        String outFile = null;
        String operation = null;

        for (String arg : args) {
            String[] split = arg.split("=");
            String key = split[0].toLowerCase();
            String value = split.length > 1 ? split[1] : null;

            if (key.equals("in_a")) {
                fileA = value;
            } else if (key.equals("in_b")) {
                fileB = value;
            } else if (key.equals("out")) {
                outFile = value;
            } else if (key.equals("op")) {
                operation = value;
            }
        }

        if (fileA == null || fileB == null || outFile == null || operation == null) {
            System.err.println("Error: Parameters 'in_a', 'in_b', 'out', and 'op' are all required.");
            System.exit(1);
        }

        try {
            // --- Core Logic ---
            System.err.println("Loading file A: " + fileA);
            Map<GeneQuad, GffLine> mapA = loadGffToMap(fileA);
            System.err.println("Loaded " + mapA.size() + " unique features.");

            System.err.println("Loading file B: " + fileB);
            Map<GeneQuad, GffLine> mapB = loadGffToMap(fileB);
            System.err.println("Loaded " + mapB.size() + " unique features.");

            Set<GeneQuad> keysA = new HashSet<>(mapA.keySet());
            Set<GeneQuad> keysB = mapB.keySet();
            
            Set<GeneQuad> resultKeys = null;

            // --- Perform Set Operation ---
            System.err.println("Performing operation: " + operation);
            if ("intersect".equalsIgnoreCase(operation)) {
                keysA.retainAll(keysB); 
                resultKeys = keysA;
            } else if ("union".equalsIgnoreCase(operation)) {
                keysA.addAll(keysB); 
                resultKeys = keysA;
            } else if ("subtract".equalsIgnoreCase(operation)) {
                keysA.removeAll(keysB); 
                resultKeys = keysA;
            } else {
                System.err.println("Error: Unknown operation '" + operation + "'. Use intersect, union, or subtract.");
                System.exit(1);
            }

            // --- Write Output ---
            System.err.println("Writing " + resultKeys.size() + " features to " + outFile);
            Map<GeneQuad, GffLine> unionMap = "union".equalsIgnoreCase(operation) ? new HashMap<>(mapA) : null;
            if (unionMap != null) { unionMap.putAll(mapB); }

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
                writer.write("##gff-version 3\n");
                writer.write("##Operation: " + operation + " on " + fileA + " and " + fileB + "\n");

                for (GeneQuad key : resultKeys) {
                    GffLine lineToWrite = null;
                    if ("union".equalsIgnoreCase(operation)) {
                        lineToWrite = unionMap.get(key);
                    } else {
                        lineToWrite = mapA.get(key);
                    }
                    writer.write(lineToWrite.toString());
                    writer.newLine();
                }
            }

            System.err.println("Done.");

        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static Map<GeneQuad, GffLine> loadGffToMap(String filePath) throws IOException {
        Map<GeneQuad, GffLine> map = new HashMap<>();
        try (BufferedReader reader = openGffReader(filePath)) {
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

    private static BufferedReader openGffReader(String filePath) throws IOException {
        if (filePath.endsWith(".gz")) {
            return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filePath))));
        }
        return new BufferedReader(new FileReader(filePath));
    }
}