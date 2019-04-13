#!/bin/bash

# Process input arguments: (1) input file of variants, (2) path to a tmp directory to store results
input_file="$1"
tmp="$2"
prefix_out=$(basename "$input_file" | cut -d '.' -f1)

# Set up header for output file
echo -e "chromosome\thg19_variant_position\treference\tvariant\tnearest_ss_dist\texon_length\texon" > "$tmp/$prefix_out.structure"

tail -n +2 "$input_file" | awk '{
                pos = $3;
                start = $6;
                stop = $7
                
                exon_length = stop - start;

                # variant is located within exon
                if(pos >= start && pos <= stop) {
                        exon = 1;
                        dist1 = pos - start;
                        dist2 = stop - pos;
                        if(dist1 <= dist2) {
                                ss_dist = dist1;
                        }
                        else {
                                ss_dist = dist2;
                        }
                }
                else {
			exon = 0;
                        if(pos > stop){
                                ss_dist = pos - stop;
                        }
                        else{
                                ss_dist = start - pos;
                        }
                }
                printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$4,$5,ss_dist,exon_length,exon);
}' | sort -k1,1 -k2,2n >> "$tmp/$prefix_out.structure"
