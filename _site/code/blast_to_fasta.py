import sys
filename = sys.argv[-1]
with open(filename) as fasta:
	for lines in fasta:
		if lines[:18]=="  <BlastOutput_db>":
			print (">"+lines[18:lines.index('.')])
		if lines[:16]=="      <Hsp_hseq>":
			new_line = (lines[16:lines.index('/')])
			print (new_line[:-1])
		else:
			continue
