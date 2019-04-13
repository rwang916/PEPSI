#!/bin/bash

score_motif()
{
        # This function calculates the change presence of putative splicing regulatory elements

        prefix=$(basename "$1")

        while read line; do
		variant_id=$(echo "$line" | cut -f1)
                chr=$(echo "$line" | cut -f2)
                pos=$(echo "$line" | cut -f3)
                ref=$(echo "$line" | cut -f4)
                alt=$(echo "$line" | cut -f5)
		
		ref_lb=$(grep "$variant_id" "$seq_file" | cut -f2)
                ref_ub=$(grep "$variant_id" "$seq_file" | cut -f3)
                refseq=$(grep "$variant_id" "$seq_file" | cut -f4)
                alt_lb=$(grep "$variant_id" "$seq_file" | cut -f5)
                alt_ub=$(grep "$variant_id" "$seq_file" | cut -f6)
                altseq=$(grep "$variant_id" "$seq_file" | cut -f7)

                printf "%s\t%s\t%s\t%s" "$chr" "$pos" "$ref" "$alt" >> "$tmp/$prefix.motif"

                for group in "ESEseq" "ESSseq"; do
                        likelihood_score=$(python "vs_motif.py" "$refseq" "$altseq" "$data/$group/$group.list" "$prefix" "$tmp" "$ref_lb" "$ref_ub" "$alt_lb" "$alt_ub" "secondary-structure")
                        printf "\t%s" "$likelihood_score" >> "$tmp/$prefix.motif"
                done

                echo "" >> "$tmp/$prefix.motif"

        done < "$1"
}

export data
export tmp
export seq_file
export -f score_motif

# Process input arguments: (1) input file of variants, (2) path to a tmp directory to store results, (3) number of threads, (4) file of sequences
input_file="$1"
tmp="$2"
threads="$3"
seq_file="$4"

# Set up working directories
data="../data"
prefix_out=$(basename "$input_file" | cut -d '.' -f1)

num_lines=$(tail -n +2 "$input_file" | wc -l)
lines_per_file=$(awk 'BEGIN {printf("%.0f",('"$num_lines"'+'"$threads"'-1)/'"$threads"')}')
tail -n +2 "$input_file" | split --lines="$lines_per_file" - "$tmp/vs_input."

# Split input file of variants for parallel processing
find "$tmp" -name "vs_input.*" -print0 | parallel -0 -n 1 -P "$threads" -I '{}' score_motif '{}'

# Combine annotations into one file
find "$tmp" -name "vs_input.*.motif" | while read file_motif; do
        cat "$file_motif" >> "$tmp/vs_input.motif.unsorted"
done

# Set up header for output file
printf "%s\t%s\t%s\t%s" "chromosome" "hg19_variant_position" "reference" "variant" > "$tmp/$prefix_out.motif"

for group in "ESEseq" "ESSseq"; do
        printf "\t%s" "$group" >> "$tmp/$prefix_out.motif"
done

echo "" >> "$tmp/$prefix_out.motif"
sort -k1,1 -k2,2n "$tmp/vs_input.motif.unsorted" >> "$tmp/$prefix_out.motif"

rm $tmp/vs_input.*
rm *_lunp
rm *.ps
