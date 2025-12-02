package prok;

import java.util.ArrayDeque;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import shared.LineParser1;
import structures.ByteBuilder;

/**
 * Converts a GFF file with custom VECTOR attributes into a TSV file suitable
 * for machine learning. Uses ByteFile/ByteStreamWriter so both ends stream and
 * gzipped inputs/outputs are handled automatically. The VECTOR attribute often
 * encodes the contig ID (which may contain commas or digits) followed by the
 * numerical feature payload, so this implementation peels numeric tokens from
 * the tail of the attribute until all feature (and optional label) columns are
 * captured.
 */
public class GffToTsv {

    private static final int FEATURE_COLUMNS = 356;
    private static final int LABEL_COLUMNS = 1; // present when true genes were available
    private static final int MAX_NUMERIC_COLUMNS = FEATURE_COLUMNS + LABEL_COLUMNS;

    public static void main(String[] args) {
        String inFile = null;
        String outFile = null;

        for (String arg : args) {
            String[] split = arg.split("=");
            String key = split[0].toLowerCase();
            String value = split.length > 1 ? split[1] : null;

            if (key.equals("in")) { inFile = value; }
            else if (key.equals("out")) { outFile = value; }
        }

        if (inFile == null || outFile == null) {
            System.err.println("Usage: java prok.GffToTsv in=<file.gff> out=<file.tsv>");
            System.exit(1);
        }

        System.err.println("Converting " + inFile + " to " + outFile);

        ByteFile reader = null;
        ByteStreamWriter writer = null;
        LineParser1 parser = new LineParser1('\t');
        ByteBuilder bb = new ByteBuilder(4096);
        long written = 0;

        try {
            reader = ByteFile.makeByteFile(inFile, true);
            writer = new ByteStreamWriter(outFile, true, false, true);
            writer.start();

            for (byte[] line = reader.nextLine(); line != null; line = reader.nextLine()) {
                if (line.length == 0 || line[0] == '#') { continue; }
                parser.set(line);
                if (parser.terms() < 9) { continue; }

                String attributes = parser.parseString(8);
                int vectorIndex = attributes.indexOf("VECTOR=");
                if (vectorIndex < 0) { continue; }

                String vectorData = attributes.substring(vectorIndex + 7);
                int semi = vectorData.indexOf(';');
                if (semi >= 0) { vectorData = vectorData.substring(0, semi); }

                ArrayDeque<String> numericTokens = extractNumericTokens(vectorData);
                if (numericTokens.size() < FEATURE_COLUMNS) {
                    System.err.println("WARNING: Skipping line due to insufficient numeric columns (found "
                            + numericTokens.size() + ")");
                    continue;
                }

                bb.clear();
                boolean first = true;
                for (String token : numericTokens) {
                    if (!first) { bb.append('\t'); }
                    bb.append(token);
                    first = false;
                }
                writer.println(bb);
                written++;
            }

            System.err.println("Conversion complete. Wrote " + written + " rows.");

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (reader != null) { reader.close(); }
            if (writer != null) { writer.poisonAndWait(); }
        }
    }

    private static ArrayDeque<String> extractNumericTokens(String vectorData) {
        ArrayDeque<String> tokens = new ArrayDeque<>(MAX_NUMERIC_COLUMNS);
        int end = vectorData.length();
        int numericCount = 0;

        while (end > 0) {
            int start = vectorData.lastIndexOf(',', end - 1);
            int tokenStart = (start < 0) ? 0 : start + 1;
            String token = vectorData.substring(tokenStart, end).trim();

            if (!token.isEmpty()) {
                try {
                    Double.parseDouble(token);
                    tokens.addFirst(token);
                    numericCount++;
                    if (numericCount == MAX_NUMERIC_COLUMNS) {
                        break;
                    }
                } catch (NumberFormatException e) {
                    break; // reached textual contig segment
                }
            }

            if (start < 0) { break; }
            end = start;
        }

        return tokens;
    }
}
