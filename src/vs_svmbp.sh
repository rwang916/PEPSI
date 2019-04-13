#!/bin/bash

svm_bp()
{
        # This is a function that calculates the change in branch point sequence strength
        # due to a mutation (SNV/small indel) using SVM-BP

        prefix=$(basename "$1")

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

		threshold=$(bc -l <<< "$ref_lb-1")
                refseq=$(echo "$wt_seq" | cut -c1-"$threshold")
                echo -e ">reference\n$refseq" > "$tmp/$prefix.ref.fa"
                python "$data/SVMBP/svm_bpfinder.py" -i "$tmp/$prefix.ref.fa" -s "Hsap" -p "$prefix" -o "$tmp" > "$tmp/$prefix.ref.svmbp"

                # If SVM-BP identifies potential branch point sequences, output the best one
                if [ -s "$tmp/$prefix.ref.svmbp" ]; then
                        svm_ref=$(perl "$data/SVMBP/calculate_best_BP_per_intron.pl" < "$tmp/$prefix.ref.svmbp" | cut -f10)
                        if [ -z "$svm_ref" ]; then
                                svm_ref="0"
                        fi
                else
                        svm_ref="0"
                fi

                threshold=$(bc -l <<< "$alt_lb-1")
                altseq=$(echo "$mut_seq" | cut -c1-"$threshold")
                echo -e ">variant\n$altseq" > "$tmp/$prefix.alt.fa"
		python "$data/SVMBP/svm_bpfinder.py" -i "$tmp/$prefix.alt.fa" -s "Hsap" -p "$prefix" -o "$tmp" > "$tmp/$prefix.alt.svmbp"

                # If SVM-BP identifies potential branch point sequences, output the best one
                if [ -s "$tmp/$prefix.alt.svmbp" ]; then
                        svm_alt=$(perl "$data/SVMBP/calculate_best_BP_per_intron.pl" < "$tmp/$prefix.alt.svmbp" | cut -f10)
                        if [ -z "$svm_alt" ]; then
                                svm_alt="0"
                        fi
                else
                        svm_alt="0"
                fi

                # Compute the difference in branch point sequence strengths between the mutated sequence and reference sequences
                svm_score=$(awk -v "svm_alt=$svm_alt" -v "svm_ref=$svm_ref" 'BEGIN{print svm_alt-svm_ref}')

                rm "$tmp/$prefix.ref.svmbp"
                rm "$tmp/$prefix.ref.fa"
                rm "$tmp/$prefix.alt.svmbp"
                rm "$tmp/$prefix.alt.fa" 
		
                echo -e "$chr\t$pos\t$ref\t$alt\t$svm_score" >> "$tmp/$prefix.svmbp"

        done < "$1"
}

export tmp
export seq_file
export data
export -f svm_bp

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
find "$tmp" -name "vs_input.*" -print0 | parallel -0 -n 1 -P "$threads" -I '{}' svm_bp '{}'

# Combine annotations into one file
find "$tmp" -name "vs_input.*.svmbp" | while read file_svmbp; do
        cat "$file_svmbp" >> "$tmp/vs_input.svmbp.unsorted"
done

# Set up header for output file
echo -e "chromosome\thg19_variant_position\treference\tvariant\tSVM_BP" > "$tmp/$prefix_out.svmbp"
sort -k1,1 -k2,2n "$tmp/vs_input.svmbp.unsorted" >> "$tmp/$prefix_out.svmbp"

rm $tmp/vs_input.*
