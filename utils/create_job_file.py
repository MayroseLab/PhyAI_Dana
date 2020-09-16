import os, sys
from subprocess import call
__author__ = 'Nomi'

sys.path.append(os.path.dirname(sys.path[0]))

from defs import *
#from utils import change_path_permissions_to_777

# qsub command and arguments
QSUB_ARGS = \
'''
#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q {queue}
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N {job_name}
#PBS -e {dirpath}
#PBS -o {dirpath}
cd $PBS_O_WORKDIR

module load python/python-anaconda3.7-itaym
{commands}
'''

#module load python/anaconda_python-3.5

QUEUE = "itaym"    #r@power9"
SH_FILE = "job.sh"
DIRECTORY = "run_{}"
PERMISSION_CMD = "chmod -R 777 ."

CMD = ''

def run_job(cmd, dirpath, sh_file=SH_FILE, job_name="job", priority=-1, chmod=False):
	error_files_path = dirpath + "error_files/"
	if not os.path.exists(error_files_path):
		os.mkdir(error_files_path)

	#if change permission
	if chmod:
		cmd = "{}\n{}".format(cmd,PERMISSION_CMD)

	#format qsub arguments
	qsub_args = QSUB_ARGS.format(queue=QUEUE,
								 commands=cmd,
								 priority=priority,
								 job_name=job_name,
								 dirpath=error_files_path)

	# create the arguments file for qsub
	with open(error_files_path + sh_file, 'w') as f:
		f.write(qsub_args)

	#call qsub command
	cmd = "qsub -p {priority} {sh_file}".format(priority=priority, sh_file=error_files_path + sh_file)
	call(cmd.split(" "))


def main(command, dirpath, sh_file, multiply_jobs, priority, job_name):
	if not multiply_jobs:
		run_job(command, dirpath, sh_file, job_name, priority)
	else:
		for i in range(1, multiply_jobs + 1):
			dir = DIRECTORY.format(i)
			os.makedirs(dir)
			os.chdir(dir)
			run_job(command, dirpath, sh_file, job_name, priority)
			os.chdir("../")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--command', '-cmd', required=False,
						default=CMD, help='command to run')
	parser.add_argument('-n', default=0, type=int,
						help='number of times to run command')
	parser.add_argument('--dirpath', '-dp', default=DIRPATH,
						required=False, help='path for sh an ER files')
	parser.add_argument('--sh_file', '-sh', default=SH_FILE,
						required=False, help='name of .sh file')
	parser.add_argument('--priority', '-p', default=-1, type=int,
						required=False, help='priority of job')
	parser.add_argument('-job_name', default="job",
						required=False, help='name of job')

	args = parser.parse_args()

	main(args.command, args.dirpath, args.sh_file, args.n, args.priority, args.job_name)

