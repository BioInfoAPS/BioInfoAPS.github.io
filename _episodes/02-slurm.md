---
title: "SLURM workload manager"
teaching: 10
exercises: 15
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# SLURM workload manager

Unlike a personal computer, a computation cluster is used by large number 
of users to run resource intensive jobs. 
The sum of the resources requested by 
those jobs often exceed the resources available in cluster. 
To deal with this issue, many high performance clusters use schedulers 
to manage resources allocated to those jobs.

UF HiperGator uses a popular scheduler called Slurm workload manager.  
Portable Batch System (PBS), Platform LSF, Moab, LoadLeveler etc. 
are examples of other schedulers.

Slurm is has three major functions:

1. Allocate exclusive and/or non-exlcudisve access to resources for a duration of time so they can perform a work,
2. Provide a framework for starting, executing, and monitoring work on allocated nodes,
3. Arbitrate contention for resources by managing a queue of pending work.

Advantage of using (SLURM) scheduler:
1. Once resources are allocated to a job, they are not taken away until job exits. If someone runs intensive
job at the same time, the resources available to the job is not changed. This improves speed and reliability
of job completion.
2. Unlike interactive window, SLURM jobs do not stop when user is not logged in.

### SLURM scripts

In order to submit a job, user must provide a scripts that will be used to specify user, 
account, time limit, memory, job output and other information.
A common way to do this is to create a submission script. 
A submission script, like other bash scripts, 
contains all the commands to be run by the job. 
However, it alsocontains extra comments which can be read by SLURM 
to determine resource requests.

To get started with SLURM, lets create new directory in your work directory and 
copy a submission script template from share/scripts directory.

```sh
$ cd /blue/general_workshop/<username>

$ mkdir slurm

$ cd slurm

$ cp ../../share/scripts/slurm_template.sh ./slurm.sh
```

Note: `.sh` is commonly used extension for shell scripts. Using a extension is not mandatory.

### Adding information to SLURM script

We have to modify some information in the template to make the provide more information
to SLURM about the job.

We can use a small text editor program called `nano` for writing to a file.

```sh
$ nano slurm.sh
```

A basic text editor will open.
```
-----------------------------------------------------------------------------------------------
 GNU nano 3.3 beta 02                      New Buffer
-----------------------------------------------------------------------------------------------
#!/bin/bash
#
#SBATCH --job-name serial_job_test    # Job name
#SBATCH --account general_workshop    # Account to run the computational task
#SBATCH --qos general_workshop        # Account allocation
#SBATCH --mail-type END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user <email_address>   # Where to send mail  
#SBATCH --ntasks 1                    # Run on a single CPU
#SBATCH --mem 10mb                    # Job memory request
#SBATCH --time 00:05:00               # Time limit hrs:min:sec
#SBATCH --output serial_test_%j.log   # Standard output and error log

-----------------------------------------------------------------------------------------------
^G Get Help     ^O WriteOut     ^R Read File     ^Y Prev Page     ^K Cut Text       ^C Cur Pos
^X Exit         ^J Justify      ^W Where Is      ^V Next Page     ^U UnCut Text     ^T To Spell
-----------------------------------------------------------------------------------------------
```

The comments beginning with `#SBATCH` tell SLURM various information about the job.

Change the &lt;email_address&gt; to your temporary UF email adress.

Lets add more lines under the last line to run a demo job.

``` 
# Return current date and time every 30 seconds for 6 times
for i in {0..5}; do printf '%s %s\n' "$(date)"; sleep 10s; done
```

Press <kbd>Ctrl</kbd>+<kbd>x</kbd> (<kbd>Cmd</kbd>+<kbd>x</kbd> in MacOS) to return to bash prompt.
Press <kbd>Y</kbd> to save the changes made to the file.

## Running a job in SLURM

### Submitting a SLURM job

To submit the job to SLURM, `sbatch` command is used.

```sh
$ sbatch slurm.sh 
Submitted batch job <jobid>
```

### Checking status of s SLURM job

You can check the status of the job using the command `squeue`. 
`-u` argument accepts &lt;username&gt; and displays all jobs by the user.
`-A` argument accepts account name and displays all jobs 
using the resources allocated to that account.

```sh
$ squeue -u <username> 
    JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
  <jobid> hpg2-comp serial_j   <user>  R       0:07      1 c29a-s2


$ squeue -A general_workshop
    JOBID PARTITION     NAME      USER  ST       TIME  NODES NODELIST(REASON)
  <jobid> hpg2-comp serial_j   <user1>  PD       0:00      1 c29a-s2
  <jobid> hpg2-comp serial_j   <user2>   R       1:12      1 c15a-s1
  <jobid> hpg2-comp serial_j   <user3>   R       0:46      1 c09a-s4
```

Under status `ST`, `R` stands for Running and `PD` stands for pending.


SLURM script: (/blue/general_workshop/share/scripts/slurm_template.sh)
```
#!/bin/bash
#SBATCH --job-name=serial_job_test    # Job name
#SBATCH --account=general_workshop    # Account to run the computational task
#SBATCH --qos=general_workshop        # Account allocation
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email_address>   # Where to send mail  
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
```