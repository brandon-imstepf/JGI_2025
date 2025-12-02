#!/bin/bash

# Parse arguments
for arg in "$@"; do
    case $arg in
        in=*)
            INPUT="${arg#*=}"
            ;;
        out=*)
            OUTPUT="${arg#*=}"
            ;;
        *)
            echo "Unknown argument: $arg"
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Usage: $0 in=<input.tsv> out=<output_clean.tsv>"
    exit 1
fi

# Remove empty lines for processing
tmp_input=$(mktemp)
# Instead of filtering out empty lines, just copy the input file to the temp file
cp "$INPUT" "$tmp_input"

# Find the most common number of columns from the first 10000 lines (excluding empty lines)
echo "Column counts in first 20 lines for debugging:"
head -n 20 "$tmp_input" | awk -F'\t' '{print NR ": " NF " columns"}'
most_common_cols=$(head -n 10000 "$tmp_input" | awk -F'\t' 'NF>0{print NF}' | sort | uniq -c | sort -nr | head -1 | awk '{print $2}')
echo "All unique column counts and their frequencies in first 10000 lines:"
head -n 10000 "$tmp_input" | awk -F'\t' 'NF>0{print NF}' | sort | uniq -c | sort -nr
# Debug: output number of rows in the entire file:
total_rows=$(wc -l < "$tmp_input")
echo "Total rows in input file: $total_rows"

#most_common_cols=354 # Just make it this for now, since we know the data is 354 columns

if [[ -z "$most_common_cols" ]]; then
    echo "ERROR: Could not determine the most common column count. Is the file empty or malformed?"
    rm "$tmp_input"
    exit 1
fi

echo "Most common column count (from first 10k lines): $most_common_cols"

# Normalize line endings and remove trailing tabs before filtering
tmp_cleaned=$(mktemp)
sed 's/\r$//' "$tmp_input" | sed 's/\t*$//' > "$tmp_cleaned"
awk -F'\t' -v cols="$most_common_cols" 'NF==cols' "$tmp_cleaned" > "$OUTPUT"
rm "$tmp_cleaned"

total_lines=$(grep -c . "$tmp_input")
good_lines=$(grep -c . "$OUTPUT")
bad_lines=$((total_lines - good_lines))

if [[ $total_lines -gt 0 ]]; then
    percent=$(awk "BEGIN {printf \"%.2f\", ($bad_lines/$total_lines)*100}")
else
    percent="0.00"
fi

echo "Removed $bad_lines corrupted lines out of $total_lines ($percent%)"
echo "Cleaned file written to $OUTPUT"

rm "$tmp_input"