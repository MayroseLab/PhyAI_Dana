import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")

from defs import *


def check_num_jobs(filepath):
	os.system("qstat -u danaazouri | wc -l > {}".format(filepath))
	time.sleep(5)
	with open(filepath, 'r') as fp:
		num_lines = fp.read()
		if int(num_lines)-5 < 300:
			return True  # can submit more
	return False



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='run SPR dataset by dataset')
	parser.add_argument('--runover', '-r', default=False, action='store_true')
	args = parser.parse_args()

	df = pd.read_csv(SUMMARY_FILES_DIR + CHOSEN_DATASETS_FILENAME)
	filepath = "{}jobs_counter.txt".format(DIRPATH)
	for i in range(len(df)):
		job_name = "SPR_for_ds_{}".format(i)
		cmd = "python " + CODE_PATH + "execute_programs/SPR_move.py " + "-dsi " + str(i)
		if args.runover:
			cmd += " -r "
		create_job_file.main(command=cmd, dirpath=DIRPATH, sh_file=job_name + ".sh", multiply_jobs=False,
							 priority=0, job_name=job_name)

		while not check_num_jobs(filepath):		# queue is full
			print("sleeping 2 mins")
			time.sleep(120)

		if i ==1 :
			exit()