---
title: "Parsing BLAST output"
teaching: 20
exercises: 15
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# File parsing

In a bioinformatics project, it is common to use several programs in series
i.e. output of a program is used as input to another program. 
There are many standard data formats such as FASTA and BAM which are 
commonly accepted by most programs. 
This makes transition from one program to another program easier.
Unfortunately, sometimes programs output non-standard or 
proprietary formats, causing disconnect in using downstream programs. 
This is where file parsing is important. 
File parsing simply means reading the file and facilitating the data 
transformation for subsequent analysis.

Note: While file parsing does not sound important, a bioinformatician has 
to spend significant time corectig file formats.

## Running BLASTn for phytogenetic tree

In next lesson, we will phylogenetically compare ‘avrBs2’ gene 
from various strains of Xanthomonas species.
The ouline of activity will be as follows:

BLAST &rarr; parse to fasta format &rarr; align sequence &rarr; phylogenetic analysis
 
To get started, we will run BLASTn in SLURM just like last section, 
however, we will work with a larger set (20) of *Xanthomonas* genomes.
You can find the genomes in the directory 
`/blue/general_workshop/share/phylogeny`. 

```sh
$ cd /blue/general_workshop/<username>

$ cp -r ../share/phylogeny ./

$ cd phylogeny

$ ls
avrBs2.fas     GEV1026.fasta     GEV1001.fasta     GEV1044.fasta     GEV1054.fasta
...
...
```

The SLURM submission script is located in `/blue/general_workshop/share/scripts/slurm_blast_xml.sh`

Warning: This is not the same script as last section. Use this script, not the last one.

```sh
$ tail -n8 ../share/scripts/slurm_blast_xml.sh
for genome in `ls *.fasta | sed 's/.fasta//g'`
do
  # Make database
  makeblastdb -in "$genome.fasta" --dbtype nucl -out "$genome"

  # Run blastn on that database
  blastn -query avrBs2.fas -db "$genome" -out $genome"_avrBs2.out" -outfmt 5 -evalue 0.001
done
```

Note the new argument `-outfmt 5`. This will output the result in XML format.


## Exercise

Submit the SLURM script. 
Do not forget to copy it to your current (`phylogeny`) directory and 
to replace your own email address first.

---

The folder should now have 20 output files with ’_avrBs2.out’ from each genome. 

```sh
$ ls *_avrBs2.out | wc -l
20
```
---

Optional: Don’t forget to check one of the output file using a text editor 
such as ‘BBedit’ for mac users and ‘Sublime text’ for windows.
The output file in this output format is full of information about the 
query and database sequence. 
Optional because Filezilla may take time!!

## BLAST XML output Vs FASTA format

Now, we can convert the BLAST `XML` output in a more common format
that can be read by multiple software to run analyses. 
Specifically, we will convert it to `FASTA` format, which we
will subsequently use to generate a phylogenetic tree.

BLAST XML format ncludes a lot of information on query and database name, 
sequence identity, statistics, alignments and so on.
Below, you can see what BLAST XML output looks like:

Note: The example below has been abridged.

```
$ cat 
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.10.1+</BlastOutput_version>
  <𝗕𝗹𝗮𝘀𝘁𝗢𝘂𝘁𝗽𝘂𝘁_𝗱𝗯>𝗫𝗲𝘂</𝗕𝗹𝗮𝘀𝘁𝗢𝘂𝘁𝗽𝘂𝘁_𝗱𝗯>
</BlastOutput>
...
...
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_def>NZ_CM002866.1 Xanthomonas euvesicatoria pv. allii CFBP 6369 chromosome, whole genome shotgun sequence</Hit_def>
  <Hit_hsps>
    <Hsp>
      <Hsp_bit-score>4170.86</Hsp_bit-score>
      <Hsp_score>2258</Hsp_score>
      <Hsp_evalue>0</Hsp_evalue>
      <Hsp_identity>2288</Hsp_identity>
      <Hsp_positive>2288</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>2303</Hsp_align-len>
      <Hsp_qseq>AAGTGCTGGCAACGCGTCCAAACACGAGCAGGCCAGGCAGACCGAGACGGATTGA</Hsp_qseq>
      <𝗛𝘀𝗽_𝗵𝘀𝗲𝗾>𝗔𝗔𝗚𝗧𝗚𝗖𝗧𝗚𝗚𝗖𝗔𝗔𝗖𝗚𝗖𝗚𝗧𝗖𝗖𝗔𝗔𝗔𝗖𝗔𝗖𝗖𝗔𝗚𝗖𝗔𝗚𝗚𝗖𝗖𝗔𝗚𝗚𝗖𝗔𝗚𝗔𝗖𝗖𝗚𝗔𝗚𝗔𝗖𝗚𝗚𝗔𝗧𝗧𝗚𝗔</𝗛𝘀𝗽_𝗵𝘀𝗲𝗾>
      <Hsp_midline>||||||||||||||||||||||||| |||||||||||||||||||||||||||||</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
```
Warning: The `.out` files you generated in last section is not in XML format.

We will have to convert it into FASTA format which only needs two line from the XML:
`BlastOutput_db` and `Hsp_hseq`. The FASTA format looks like:

```
> Xeu
AAGTGCTGGCAACGCGTCCAAACACCAGCAGGCCAGGCAGACCGAGACGGATTGA
```

### Bash script to parse XML to FASTA

We will now use a bash script to parse the BLAST output to 
extract the file as a multiFASTA file.
THe script is located in path `/blue/general_workshop/share/scripts/blast2fasta.sh`

```sh
$ cp ../../share/scripts/blast2fasta.sh ./   

$ cat blast2fasta.sh
#!/bin/bash
while read line
do
        [[ `grep "BlastOutput_db" <<< $line` ]] && echo -n ">" <<< $line
        [[ `grep -E "BlastOutput_db|Hsp_hseq" <<< $line` ]] && sed -n "s:.*>\(.*\)</.*:\1:p" <<< $line
done
```

What this script does is loop line by line over the file, 
and look for keywords `BlastOutput_db` and `Hsp_hseq`.
If it find the keywords, then it will extract the contents
between those tags and output in appropriate format.

Tip: `sed` is replacing everything except necessary information with empty string.

### Runnig the bash script

We can run the script now.
We can call the script in two ways

- `./<script_name>`
- `sh <scriptname>`

```sh
$ cat *_avrBs2.out | ./blast2fasta.sh > avrBs2_all_genomes.fas

$ cat avrBs2_all_genomes.fas
>GEV1026
AAGTGCTGGCAACGCGTCCAAACAAGGCCTGCGCCGCACGCCTGCCAGCGCGCGCAACGCAGGCATCGTTTCGCATCCGGG
CGGTACTTTTCGCCTAATTTGCCAATTGTCATATGCCACGCGCTTTACTGGCCGCCCGCCGCGTTTTCGAGGTCATCATGC
GCATCGGTCCTCTGCAACCTTCTATCGCGCACACTGC
...
...
>GEV1001
GTGCTGGCAACGCGTCCAAACAAGGCCTGCGCCGCACGCCTGCCAGCGCGCGTAACGCAGGCATCGTTTCGCATCGGCGGG
CGGTACTTTTCGCCTAATTTGCCAATTGTCAGATGCCACGCGCTTTACTGGCCGCCCGACGCGTTTTCGAGGTCATCATGC
GCATCGGTCCTCTGCAACCTTCTATCGCGCA
...
...
```

Tip: We are not using SLURM for this job because it is fast and not resource intensive.

Note: The example above is for demonstration purpose and the same as the
output of the script.

Now we have a multifasta file with sequence id for each genome database 
we used along with their respective sequence. 
We can now use this fasta formatted file to run sequence alignment and 
subsequently construct a tree.