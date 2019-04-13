#!/bin/bash

# Process input arguments: (1) input file of variants, (2) path to a tmp directory to store results, (3) mode (training/test)
input_file="$1"
tmp="$2"
mode="$3"

# Set up working directories
data="../data"

# Setup new header for test file
echo -e "ID\tchromosome\thg19_variant_position\treference\tvariant\texon_start\texon_end\tHepG2_delta_logit_psi" > "$tmp/vs_""$mode""_set.tsv"

# Annotate test file with exon orientation information
tail -n +2 "$input_file" | awk '{
	ref_psi = $9/100;
	variant_psi = ($8+$9)/100;

	# Restrict reference PSI to range [1e-4, 1-1e-4]
	if(ref_psi <= 1e-4) {
		ref_psi = 1e-4;
	}

	if(ref_psi >= 1-1e-4){
		ref_psi = 1 - 1e-4;
	}

	# Restrict variant PSI to range [1e-4, 1-1e-4]
	if(variant_psi <= 1e-4) {
		variant_psi = 1e-4;
	}

	if(variant_psi >= 1-1e-4){
		variant_psi = 1 - 1e-4;
	}
                
	delta_logit_psi = log(variant_psi/(1-variant_psi)) - log(ref_psi/(1-ref_psi));

	printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,delta_logit_psi)
        
}' | sort -k2,2 -k3,3n >> "$tmp/vs_""$mode""_set.tsv"

