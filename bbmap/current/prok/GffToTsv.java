package prok;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.zip.GZIPInputStream;

/**
 * Converts a GFF file with custom VECTOR attributes into a TSV file
 * suitable for machine learning. This version intelligently removes the
 * entire text-based contig ID from the vector, even if it contains commas.
 */
public class GffToTsv {

    public static void main(String[] args) {
        String inFile = null;
        String outFile = null;
        String label = "0"; // Default label for the 'match' column

        for (String arg : args) {
            String[] split = arg.split("=");
            String key = split[0].toLowerCase();
            String value = split.length > 1 ? split[1] : null;

            if (key.equals("in")) inFile = value;
            else if (key.equals("out")) outFile = value;
            else if (key.equals("label")) label = value;
        }

        if (inFile == null || outFile == null) {
            System.err.println("Usage: java prok.GffToTsv in=<file.gff> out=<file.tsv> label=<value>");
            System.exit(1);
        }

        try (BufferedReader reader = openReader(inFile);
             BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
            
            System.err.println("Converting " + inFile + " to " + outFile + " with label " + label);

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;

                String[] parts = line.split("\t");
                if (parts.length < 9) continue;

                String attributes = parts[8];
                int vectorIndex = attributes.indexOf("VECTOR=");
                
                if (vectorIndex != -1) {
                    String fullVector = attributes.substring(vectorIndex + 7);
                    
                    String[] vectorParts = fullVector.split(",");
                    int numericalStartIndex = -1;

                    // Find the first part of the vector that is a valid number.
                    for (int i = 0; i < vectorParts.length; i++) {
                        try {
                            // We only need to check if it's a valid double.
                            Double.parseDouble(vectorParts[i].trim());
                            numericalStartIndex = i;
                            break; // Found the first number, stop searching.
                        } catch (NumberFormatException e) {
                            // This part is not a number, so it must be part of the text ID. Continue.
                        }
                    }

                    if (numericalStartIndex != -1) {
                        // Rebuild the string as a tab-separated line, starting from the first numerical part.
                        StringBuilder tsvBuilder = new StringBuilder();
                        for (int i = numericalStartIndex; i < vectorParts.length; i++) {
                            if (i > numericalStartIndex) {
                                tsvBuilder.append('\t');
                            }
                            tsvBuilder.append(vectorParts[i].trim());
                        }
                        
                        // The tsvBuilder now contains the full line, including the label from the vector.
                        // The 'label' argument from the command line is ignored as it's redundant.
                        writer.write(tsvBuilder.toString());
                        writer.newLine();
                    }
                    // If no numerical part is found in the vector, we simply skip writing the line.
                }
            }
            System.err.println("Conversion complete.");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static BufferedReader openReader(String filePath) throws IOException {
        InputStream is = new FileInputStream(filePath);
        if (filePath.endsWith(".gz")) {
            is = new GZIPInputStream(is);
        }
        return new BufferedReader(new InputStreamReader(is));
    }
}
