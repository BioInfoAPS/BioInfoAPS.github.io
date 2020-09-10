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

Unlike a personal computer, a computation cluster is used by large number 
of users to run resource intensive jobs. 
The sum of the resources requested by 
those jobs often exceed the resources available in cluster. 
To deal with this issue, many high performance clusters use schedulers 
to manage resources allocated to those jobs.

[UF Hipergator](https://www.rc.ufl.edu/services/hipergator/) uses a popular scheduler called [SLURM](https://slurm.schedmd.com/documentation.html) workload manager.  
Portable Batch System (PBS), Platform LSF, Moab, LoadLeveler etc. 
are examples of other schedulers.

SLURM is has three major functions:

1. Allocate exclusive and/or non-exlcudisve access to resources for a duration of time so they can perform a work,
2. Provide a framework for starting, executing, and monitoring work on allocated nodes,
3. Arbitrate contention for resources by managing a queue of pending work.

Advantage of using (SLURM) scheduler:
1. Once resources are allocated to a job, they are not taken away until job exits. If someone runs intensive
job at the same time, the resources available to the job is not changed. This improves speed and reliability
of job completion.
2. Unlike interactive window, SLURM jobs do not stop when user is not logged in.

### SLURM scripts

In order to submit a job, user must provide a scripts that will specify user, 
account, time limit, memory, job output, and other information.
A common way to do this is to create a submission script. 
A submission script, like other bash scripts, 
contains all the commands to be run during the job. 
However, it also contains extra comments which can be read by SLURM 
to determine resource requests.

To get started with SLURM, lets create new directory in your work directory and 
copy a submission script template from share/scripts directory.

~~~
$ cd /blue/general_workshop/<username>

$ mkdir slurm

$ cd slurm

$ cp ../../share/scripts/slurm_template.sh ./slurm.sh
~~~
{: .language-bash}

Note: `.sh` is commonly used extension for shell scripts. Using a extension is not mandatory.

### Adding information to SLURM script

We have to modify some information in the template to make the provide more information
to SLURM about the job.

We can use a small text editor program called `nano` for writing to a file.
This will open a basic text editor.

~~~
$ nano slurm.sh
~~~
{: .language-bash}

~~~
-----------------------------------------------------------------------------------------------
 GNU nano 3.3 beta 02                     File: slurm.sh
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

# Return current date and time every 3 seconds for 5 times
for i in {0..4}; do printf '%s %s\n' "$(date)"; sleep 3s; done


-----------------------------------------------------------------------------------------------
^G Get Help     ^O WriteOut     ^R Read File     ^Y Prev Page     ^K Cut Text       ^C Cur Pos
^X Exit         ^J Justify      ^W Where Is      ^V Next Page     ^U UnCut Text     ^T To Spell
-----------------------------------------------------------------------------------------------
~~~
{: .terminal}

The comments beginning with `#SBATCH` tell SLURM various information about the job.
The acutal commands to run appear after these comments. In this case, it just returns
current datetime at 3s interval.

Change the &lt;email_address&gt; to your UF email address.
Once you are done, press <kbd>Ctrl</kbd>+<kbd>x</kbd> 
(<kbd>Cmd</kbd>+<kbd>x</kbd> in MacOS) to return to bash prompt.
Press <kbd>Y</kbd> to save the changes made to the file.

> ## Editing in nano
> **nano** is a commandline editor. You can only move your cursor with 
> arrow keys: <kbd>↑</kbd>, <kbd>↓</kbd>, <kbd>←</kbd> and <kbd>→</kbd>.
> Clicking with mouse does not change the position of the cursor. 
> Be careful, you may be editing in wrong place.
>
> If you accidentally edited in wrong place, exit nano, 
> delete the script `slurm.sh` and copy again from share directory.
> Do not forget to edit the &lt;email_address>&gt;
{: .caution}

## Running a job in SLURM

### Submitting a SLURM job

To submit the job to SLURM, `sbatch` command is used.

~~~
$ sbatch slurm.sh 
~~~
{: .language-bash}

~~~
Submitted batch job <jobid>
~~~
{: .output}

### Checking status of a SLURM job

You can check the status of the job using the command `squeue`. 
`-u` argument accepts &lt;username&gt; and displays all jobs by the user.
`-A` argument accepts account name and displays all jobs 
using the resources allocated to that account.

~~~
$ squeue -u <username> 
~~~
{: .language-bash}

~~~
    JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
  <jobid> hpg2-comp serial_j   <user>  R       0:07      1 c29a-s2
~~~
{: .output}

> If you do not see your job, it may have already been completed.
> Run the job again and check within a minute.
> Also, check you email to see if you got any message from HiperGator.
{: .tips}

~~~
$ squeue -A general_workshop
~~~
{: .language-bash}

~~~
    JOBID PARTITION     NAME      USER  ST       TIME  NODES NODELIST(REASON)
  <jobid> hpg2-comp serial_j   <user1>  PD       0:00      1 c29a-s2
  <jobid> hpg2-comp serial_j   <user2>   R       1:12      1 c15a-s1
  <jobid> hpg2-comp serial_j   <user3>   R       0:46      1 c09a-s4
~~~
{: .output}

> ## Understanding Job Status
> Under status `ST`, `R` stands for Running and `PD` stands for pending.
> If the job is pending, a reason may be provided in last column. Eg:
> - **None**: Just taking a while before running.
> - **Priority**: Higher priority jobs exist in this partition.
> - **QOSMaxCpuPerUserLimit**: The user is already using max 
> number of CPU that they are allowed to use.
{: .tips}

### Checking the output

The SLURM submission script containas a line 
`#SBATCH --output serial_test_%j.log`. Thus the output for this job
with be in the file `serial_test_<jobid>.log`.

~~~
$ ls
~~~
{: .language-bash}

~~~
serial_test_<jobid>.log     slurm.sh
~~~
{: .output}

~~~
$ cat serial_test_<jobid>.log
~~~
{: .language-bash}

~~~
Tue Sep 15 02:04:05 EDT 2020
Tue Sep 15 02:04:15 EDT 2020
Tue Sep 15 02:04:25 EDT 2020
Tue Sep 15 02:04:35 EDT 2020
Tue Sep 15 02:04:45 EDT 2020
Tue Sep 15 02:04:55 EDT 2020
~~~
{: .output}

> ## Autocomplete in terminal
> To autocomplete a file or directory name, press <kbd>Tab</kbd> button. 
> Names are autocompleted until there are conflicts (e.g. files with same prefixes).
> In case of conflict, press tab two time to view list of files and folders (equivalent to `ls`).
> To try autocomplete, type `cat se` and press <kbd>Tab</kbd>.
{: .tips}
