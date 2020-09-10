#!/bin/bash
#SBATCH --job-name RAxML              # Job name
#SBATCH --account general_workshop    # Account to run the computational task
#SBATCH --qos general_workshop        # Account allocation
#SBATCH --mail-type END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user <email_address>   # Where to send mail	
#SBATCH --ntasks 1                    # Run on a single CPU
#SBATCH --mem 2gb                     # Job memory request
#SBATCH --time 00:10:00               # Time limit hrs:min:sec
#SBATCH --output RAxML_%j.log   # Standard output and error log
pwd; hostname; date

# Load RAxML
module load raxml

# Make tree with RAxML
raxmlHPC -d -p 12345 -m GTRGAMMAI -s avrBs2_all_genomes_aligned.fas -n avrBs2_tree

# Success message
echo 'ML tree output in RAxML_bestTree.avrBs2_tree'
