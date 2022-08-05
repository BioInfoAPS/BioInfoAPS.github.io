#!/bin/bash
#SBATCH --job-name phylogenetic_tree  # Job name
#SBATCH --account general_workshop    # Account to run the computational task
#SBATCH --qos general_workshop        # Account allocation
#SBATCH --mail-type END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user <email_address>   # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=00:10:00               # Time limit hrs:min:sec
#SBATCH --output tree_pipeline_%j.log   # Standard output and error log
pwd; hostname; date

# Load programs
module load ncbi_blast
module load mafft
module load raxml

# Loop over each genome
for genome in `ls *.fasta | sed 's/.fasta//g'`
do
  # Make database
   makeblastdb -in "$genome.fasta" -dbtype nucl -out "$genome"

  # Run blastn on that database
   blastn -query avrBs2.fas -db "$genome" -out $genome"_avrBs2.out" -outfmt 5 -evalue 0.001
done

# Combine all outputs and parse into multiFASTA file
cat *_avrBs2.out | ./blast2fasta.sh > avrBs2_all_genomes.fas

# Align sequences with MAAFT
mafft avrBs2_all_genomes.fas > avrBs2_all_genomes_aligned.fas

# Make tree with RAxML
raxmlHPC -d -p 12345 -m GTRGAMMAI -s avrBs2_all_genomes_aligned.fas -n avrBs2_tree
