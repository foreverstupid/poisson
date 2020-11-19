#!/bin/python3

import sys
import os
import time
import subprocess as sp

ref = f'''
    Usage:
        {__file__} <proc count> <time> <x count> <y count> <max iter> <frequency>

    <proc count>:
        The count of required processes for solving.

    <time>:
        The maximum time for the required job.
    
    <x count>:
        The count of grid nodes along X-axis.

    <y count>:
        The count of grid nodes along Y-axis.

    <max iter>:
        Maximum count of iterations.

    <frequency>:
        Iteration logging frequency.

    The precision for solution is 0.05.

    NOTE: this script hardly depends on the mpisubmit.pl script, so
    if it differs this script may not work.
'''

if (len(sys.argv) != 7 or sys.argv[1] == "-h" or sys.argv[1] == "--help"):
    print(ref)
    sys.exit(0)



def execute(prog, params):
    '''Executes a program in a shell.'''

    command = prog + " " + params
    p = sp.Popen(
        command,
        stdout = sp.PIPE,
        stderr = sp.PIPE,
        shell = True,
        universal_newlines = True
    )

    out, err = p.communicate()
    exit_code = p.wait()

    return exit_code, out, err



def start_job(args):
    '''Starts a new job on HPC.'''

    command = "mpisubmit.pl"
    args = \
        f"--stdout {args.out_file} " + \
        f"--stderr {args.err_file} " + \
        f"-p {args.procs} " + \
        f"-w {args.time} " + \
        "poisson -- " + \
        f"{args.xcnt} {args.ycnt} {args.eps} " + \
        f"{args.max_iter} {args.init_iter} {args.freq} " + \
        f"{args.data_dir} "

    exit_code, out, err = execute(command, args)
    print(f"\n\nExit code of submitting: {exit_code}\nOut:\n{out}\nErr:\n{err}")
    if not "is submitted to" in out:
        print("The job is not submitted. The soving process is aborted")
        sys.exit(1)



def wait_for_job(args):
    '''Waits until the given job is finished.'''

    i = 0
    while True:
        if (os.path.isfile(args.out_file)):
            print("Job is done")
            return
        else:
            print(f"[[{i}]]: Sleeping for 300 seconds...")
            time.sleep(300)
            i += 1



def is_solved(out_dir):
    '''Checks whether the problem is fully solved.'''

    return all([
        os.path.isfile(os.path.join(out_dir, subdir, "solution.csv"))
        for subdir in os.listdir(out_dir)
    ])



def get_initial_iter(out_dir):
    '''Gets an initial iteration for the next job.'''

    iters = []
    subdirs = [
        os.path.join(out_dir, sub)
        for sub in os.listdir(out_dir)
    ]

    for file in os.listdir(subdirs[0]):
        if all([os.path.isfile(os.path.join(sub, file)) for sub in subdirs]):
            parts = file.split('_')
            if (len(parts) > 1):
                iters.append(int(parts[1].split('.')[0]))

    return max(iters)




class JobArgs:
    '''Represents arguments for a job.'''
    
    def __init__(self, args, init):
        self.procs = int(args[1])
        self.time = args[2]
        self.xcnt = int(args[3])
        self.ycnt = int(args[4])
        self.eps = 0.05
        self.max_iter = int(args[5])
        self.init_iter = init
        self.freq = int(args[6])
        
        tmp = f"{self.procs}_{self.xcnt}_{self.ycnt}"
        self.data_dir = f"out_{tmp}"
        self.out_file = f"{tmp}_{self.init_iter}.out"
        self.err_file = f"{tmp}_{self.init_iter}.err"



solved = False
init_iter = 0
attempt = 0

while (not solved):
    job_args = JobArgs(sys.argv, init_iter)
    if not os.path.isdir(job_args.data_dir):
        os.mkdir(job_args.data_dir)

    print(f"-----------------------------[{attempt}]------------------------------------")
    print(f"Init iteration: {init_iter}")
    start_job(job_args)
    wait_for_job(job_args)
    solved = is_solved(job_args.data_dir)
    if (not solved):
        init_iter = get_initial_iter(job_args.data_dir)

    attempt += 1
