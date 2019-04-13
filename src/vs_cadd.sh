#!/bin/bash

annotate_cadd()
{
        # This is a function for annotating variants based on the following CADD features:
        #       1. Percent GC in a window of +/- 75 bp
        #       2. Primate PhastCons conservation score
        #       3. Mammalian PhastCons conservation score
        #       4. Vertebrate PhastCons conservation score
        #       5. Primate PhyloP score
        #       6. Mammalian PhyloP score
        #       7. Vertebrate PhyloP score
        #       8. Neutral evolution score defined by GERP++
        #       9. Rejected Substitution' score defined by GERP++
        #       10. Background selection score
        #       11. Mutability index from Michaelson, J.J. et al. Cell 2012
        #       12. fitCons score

        prefix=$(basename "$1")
        field_string=$(cut -f2 "$feature_list" | tr '\n' ',' | sed 's/,$//')

        while read line; do
                chr=$(echo "$line" | cut -f2 | sed 's/^chr//')
                pos=$(echo "$line" | cut -f3)
                ref=$(echo "$line" | cut -f4)
                alt=$(echo "$line" | cut -f5)
                tabix "$cadd/ExAC_r0.3_inclAnno.tsv.gz" "$chr":"$pos"-"$pos" | awk -v "ref=$ref" -v "alt=$alt" '{
                        if (length(ref) == length(alt) && length(ref) > 1) {
                                ref_allele = substr(ref,1,1);
                                alt_allele = substr(alt,1,1);
                        }
                        else {
                                ref_allele = ref;
                                alt_allele = alt;
                        }
                        if ($3 == ref_allele && $5 == alt_allele) {
                                print;
                        }
                }' | head -n 1 | cut -f"$field_string" | sed -e "s/^/$chr\t$pos\t$ref\t$alt\t/" >> "$tmp/$prefix.cadd"
        done < "$1"
}

export feature_list
export cadd
export tmp
export -f annotate_cadd

# Process input arguments: (1) input file of variants, (2) list of CADD v1.3 features, (3) path to a tmp directory to store results, (4) number of threads
input_file="$1"
feature_list="$2"
tmp="$3"
threads="$4"

# Set up working directories
data="../data"
cadd="$data/cadd"

prefix_out=$(basename "$input_file" | cut -d '.' -f1)

num_lines=$(tail -n +2 "$input_file" | wc -l)
lines_per_file=$(awk 'BEGIN {printf("%.0f",('"$num_lines"'+'"$threads"'-1)/'"$threads"')}')

# Split input file of variants for parallel processing
tail -n +2 "$input_file" | split --lines="$lines_per_file" - "$tmp/vs_input."
find "$tmp" -name "vs_input.*" -print0 | parallel -0 -n 1 -P "$threads" -I '{}' annotate_cadd '{}'

# Combine all annotations into a single file
find "$tmp" -name "vs_input.*.cadd" | while read cadd; do
        cat "$cadd" >> "$tmp/vs_input.cadd.unsorted"
done

# Set up header for output file
cut -f1 "$feature_list" | tr '\n' '\t' | sed -e 's/\t$/\n/' | sed -e 's/^/chromosome\thg19_variant_position\treference\tvariant\t/' > "$tmp/$prefix_out.cadd"
sed 's/^/chr/' "$tmp/vs_input.cadd.unsorted" | sort -k1,1 -k2,2n >> "$tmp/$prefix_out.cadd"

rm $tmp/vs_input.*
