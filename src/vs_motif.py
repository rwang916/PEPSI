#!/usr/bin/python

# This is a script for calculating the impact of a variant on putative
# splicing regulatory element activity

# Input Parameters:
#       1. reference sequence
#       2. mutated sequence
#       3. list of putative regulatory motifs
#       4. file prefix
#       5. path to tmp directory

import sys
import re
import os
import subprocess

def get_score(motif, sequence, lb, ub, prefix, tmpdir, mode):
	score = 0.0
	search_string = "(?=(" + motif + "))"
	regions = [(m.start(1),m.end(1)) for m in re.finditer(search_string, sequence[lb-1:ub])]
	motif_size = len(motif)
	w_size = 80
	l_size = 40

	for pair in regions:
		lower_bound = pair[0]
		upper_bound = pair[1]
		rna_plfold_cmd = "RNAplfold -W " + str(w_size) + " -L " + str(l_size) + " -u " + str(motif_size) + " < " + tmpdir + "/" + prefix + "." + mode + ".fa"
		get_prob_cmd = "grep -P \"^" + str(upper_bound) + "\t\" " + prefix + "." + mode + "_lunp | cut -f" + str(upper_bound - lower_bound + 1)
		subprocess.call(rna_plfold_cmd, shell=True)
		prob_proc = subprocess.Popen(get_prob_cmd, stdout=subprocess.PIPE, shell=True)
		prob = float(prob_proc.stdout.read())
		score = score + prob
        
	return score

def occurrences(string, sub):
	count = start = 0
	while True:
		start = string.find(sub, start) + 1
		if start > 0:
			count+=1
		else:
			return count

# Process input parameters
refseq = sys.argv[1].strip()
altseq = sys.argv[2].strip()
profile = open(sys.argv[3])
prefix = sys.argv[4].strip()
tmpdir = sys.argv[5].strip()
ref_lb = int(sys.argv[6].strip())
ref_ub = int(sys.argv[7].strip())
alt_lb = int(sys.argv[8].strip())
alt_ub = int(sys.argv[9].strip())
config = sys.argv[10].strip()

ref_score = 0.0
alt_score = 0.0

cut_lb = max(ref_lb - 41, 0)
cut_ub = min(ref_ub + 40, len(refseq))

f = open(tmpdir + "/" + prefix + ".ref.fa", "a")
f.write(">" + prefix + ".ref\n" + refseq[cut_lb:cut_ub])
f.close()

cut_lb = max(alt_lb - 41, 0)
cut_ub = min(alt_ub + 40, len(altseq))

g = open(tmpdir + "/" + prefix + ".alt.fa", "a")
g.write(">" + prefix + ".alt\n" + altseq[cut_lb:cut_ub])
g.close()

for line in profile:
        motif = line.strip()
	if config == "secondary-structure":
		ref_score = ref_score + get_score(motif, refseq, ref_lb, ref_ub, prefix, tmpdir, "ref")
		alt_score = alt_score + get_score(motif, altseq, alt_lb, alt_ub, prefix, tmpdir, "alt")
	else:
		ref_score = ref_score + occurrences(refseq[ref_lb-1:ref_ub], motif)
		alt_score = alt_score + occurrences(altseq[alt_lb-1:alt_ub], motif)

os.remove(tmpdir + "/" + prefix + ".ref.fa")
os.remove(tmpdir + "/" + prefix + ".alt.fa")

print alt_score-ref_score
