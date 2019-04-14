#!/bin/bash

score_3ss()
{
        # This is a function that calculates the change in 3'-splice site strength
        # due to a mutation (SNV/small indel) using MaxEntScan::3ss

        prefix=$(basename "$1")
	src_dir=$(pwd)

        while read line; do
		variant_id=$(echo "$line" | cut -f1)
                chr=$(echo "$line" | cut -f2)
                pos=$(echo "$line" | cut -f3)
                ref=$(echo "$line" | cut -f4)
                alt=$(echo "$line" | cut -f5)
	
		ref_lb=$(grep "$variant_id" "$seq_file" | cut -f2)
                ref_ub=$(grep "$variant_id" "$seq_file" | cut -f3)
                wt_seq=$(grep "$variant_id" "$seq_file" | cut -f4)
                alt_lb=$(grep "$variant_id" "$seq_file" | cut -f5)
                alt_ub=$(grep "$variant_id" "$seq_file" | cut -f6)
                mut_seq=$(grep "$variant_id" "$seq_file" | cut -f7)

                # Extract relevant bounds for running MaxEntScan:3ss
                # [20 bases in intron][3 bases in exon]
		
		lower_bound=$(bc -l <<< "$ref_lb-20")
                upper_bound=$(bc -l <<< "$ref_lb+2")
                
		refseq=$(echo "$wt_seq" | cut -c"$lower_bound"-"$upper_bound")
		cd "$data/MaxEntScan"
                mes_ref=$(echo "$refseq" | perl "score3.pl" - | cut -f2)

                lower_bound=$(bc -l <<< "$alt_lb-20")
                upper_bound=$(bc -l <<< "$alt_lb+2")
                altseq=$(echo "$mut_seq" | cut -c"$lower_bound"-"$upper_bound")
                mes_alt=$(echo "$altseq" | perl "score3.pl" - | cut -f2)

                mes_score=$(awk -v "mes_alt=$mes_alt" -v "mes_ref=$mes_ref" 'BEGIN{print mes_alt-mes_ref}')
		cd "$src_dir"
		
                echo -e "$chr\t$pos\t$ref\t$alt\t$mes_score" >> "$tmp/$prefix.3ss"

        done < "$1"
}

export tmp
export seq_file
export data
export -f score_3ss

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
find "$tmp" -name "vs_input.*" -print0 | parallel -0 -n 1 -P "$threads" -I '{}' score_3ss '{}'

# Combine annotations into one file
find "$tmp" -name "vs_input.*.3ss" | while read file_3ss; do
        cat "$file_3ss" >> "$tmp/vs_input.3ss.unsorted"
done

# Set up header for output file
echo -e "chromosome\thg19_variant_position\treference\tvariant\tMaxEntScan_3ss" > "$tmp/$prefix_out.3ss"
sort -k1,1 -k2,2n "$tmp/vs_input.3ss.unsorted" >> "$tmp/$prefix_out.3ss"

rm $tmp/vs_input.*
