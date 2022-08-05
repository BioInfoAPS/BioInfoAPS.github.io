#!/bin/bash
#SBATCH --job-name serial_job_test    # Job name
#SBATCH --account general_workshop    # Account to run the computational task
#SBATCH --qos general_workshop        # Account allocation
#SBATCH --mail-type END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user <email_address>   # Where to send mail  
#SBATCH --ntasks 1                    # Run on a single CPU
#SBATCH --mem 10mb                    # Job memory request
#SBATCH --time 00:05:00               # Time limit hrs:min:sec
#SBATCH --output serial_test_%j.log   # Standard output and error log

# Return current date and time every 30 seconds for 6 times
for i in {0..5}; do printf '%s %s\n' "$(date)"; sleep 10s; done