#!/bin/bash

annotate_variants()
{
        # This is a function for annotating variants using the following features:
        #       1. CADD features (conservation, GC content)
        #       2. exon length, distance to nearest splice site
        #       3. disruption of 5'-splice site (MaxEntScan::5ss)
        #       4. disruption of 3'-splice site (MaxEntScan::3ss)
        #       5. disruption of branch point sequence (SVM-BP)
        #       6. disruption of splicing regulatory elements

        mode="$2"
	
        # Process input file by including logit transformed HepG2 delta PSI values
        bash "vs_process_input.sh" "$1" "$outdir/tmp" "$mode"

        # Annotate variants using features from CADD version 1.3
        bash "vs_cadd.sh" "$outdir/tmp/vs_""$mode""_set.tsv" "$data/cadd/vs_cadd-1.3_feature.list" "$outdir/tmp" "$threads"

        # Annotate variants by transcript structure features (exon length, distance to nearest splice splice site, etc.)
        bash "vs_structure.sh" "$outdir/tmp/vs_""$mode""_set.tsv" "$outdir/tmp"

        # Annotate variants using MaxEntScan::5ss
        bash "vs_5ss.sh" "$outdir/tmp/vs_""$mode""_set.tsv" "$outdir/tmp" "$threads" "$3"

        # Annotate variants using MaxEntScan::3ss
        bash "vs_3ss.sh" "$outdir/tmp/vs_""$mode""_set.tsv" "$outdir/tmp" "$threads" "$3"

        # Annotate variants using SVMBP
        bash "vs_svmbp.sh" "$outdir/tmp/vs_""$mode""_set.tsv" "$outdir/tmp" "$threads" "$3"

        # Annotate variants based on gain/loss of splicing regulatory elements weighed by probability of single-strandedness
        bash "vs_motif.sh" "$outdir/tmp/vs_""$mode""_set.tsv" "$outdir/tmp" "$threads" "$3"

        # Combine annotations into one large table
        cat "$outdir/tmp/vs_""$mode""_set.tsv" > "$outdir/tmp/vs_""$mode""_set.table"

        for group in "cadd" "structure" "5ss" "3ss" "svmbp" "motif"; do
                cut -f5- "$outdir/tmp/vs_""$mode""_set.$group" | paste "$outdir/tmp/vs_""$mode""_set.table" - > "$outdir/tmp/vs_""$mode""_set.table.tmp"
                rm "$outdir/tmp/vs_""$mode""_set.table"
                mv "$outdir/tmp/vs_""$mode""_set.table.tmp" "$outdir/tmp/vs_""$mode""_set.table"
        done

        cut -f1,8- "$outdir/tmp/vs_""$mode""_set.table" | sed 's/\r//g' > "$outdir/vs_""$mode""_set.table"

}

export outdir
export data
export threads
export -f annotate_variants

# Process input arguments: (1) training file, (2) test file, (3) number of threads, (4) file of training set sequences, (5) file of test set sequences
training_file="$1"
test_file="$2"
threads="$3"
training_seq_file="$4"
test_seq_file="$5"

# Set up working directories
data="../data"
testdir="../test"
mkdir -p "$testdir"

# Construct output directory
runtime=$(date +"%Y%m%d%H%M%S%N")
outdir="../test/vs-$runtime"
mkdir -p "$outdir"
mkdir -p "$outdir/tmp"

# Annotate variants in training and test files:
annotate_variants "$training_file" "training" "$training_seq_file"
annotate_variants "$test_file" "test" "$test_seq_file"

# Train random forest model with the following input parameters:
#       (1) training set of annotated variants
#       (2) test set of annotated variants
#       (3) path to file storing model predictions
#       (4) path to file storing features ranked by importance

Rscript "random_forest.R" "$outdir/vs_training_set.table" "$outdir/vs_test_set.table" "$outdir/vs_test_set.table.predictions" "$outdir/vs_model.table.varImp"

sed -i 's/"//g' "$outdir/vs_test_set.table.predictions"
sed -i 's/"//g' "$outdir/vs_model.table.varImp"

echo -e "Variant_ID\tExperimental_dPSI\tPredicted_dPSI" > "$outdir/vs_test_set.table.predictions.tmp"

tail -n +2 "$test_file" | awk '{
        if(FNR==NR) {
                delta_psi[$1]=$8;
                ref_psi[$1]=$9;
        }
        else {
                x=ref_psi[$1];
                experimental_delta_psi=delta_psi[$1];
                predicted_delta_psi=(x*(100-x)*(exp($2)-1))/(100-x+x*exp($2));
                printf("%s\t%s\t%s\n",$1,experimental_delta_psi,predicted_delta_psi)
        }
}' - "$outdir/vs_test_set.table.predictions" >> "$outdir/vs_test_set.table.predictions.tmp"

rm "$outdir/vs_test_set.table.predictions"
mv "$outdir/vs_test_set.table.predictions.tmp" "$outdir/vs_test_set.table.predictions"

sort -k2,2rg "$outdir/vs_model.table.varImp" > "$outdir/vs_model.table.varImp.tmp"
rm "$outdir/vs_model.table.varImp"
mv "$outdir/vs_model.table.varImp.tmp" "$outdir/vs_model.table.varImp"

