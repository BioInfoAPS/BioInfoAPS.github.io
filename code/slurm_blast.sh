#!/bin/bash
#SBATCH --job-name blast                    # Job name
#SBATCH --account general_workshop          # Account to run the computational task
#SBATCH --qos general_workshop              # Account allocation
#SBATCH --mail-type END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user <email_address>         # Where to send mail
#SBATCH --ntasks 1                          # Run on a single CPU
#SBATCH --mem 1gb                           # Job memory request
#SBATCH --time 00:05:00                     # Time limit hrs:min:sec
#SBATCH --output blast_job_%j.log           # Standard output and error log
pwd; hostname; date

# Load blast
module load ncbi_blast

# Loop over each genome
for genome in `ls *.fasta | sed 's/.fasta//g'`
do
  # Make database
  makeblastdb -in "$genome.fasta" -dbtype nucl -out "$genome"

  # Run blastn on that database
  blastn -query avrBs2.fas -db "$genome" -out $genome"_avrBs2.out" -evalue 0.001
done
