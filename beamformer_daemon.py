import os
import sys
from os.path import isdir, isfile, join
from subprocess import Popen
import subprocess
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def execute_cmd(cmd):
	try:
		proc = Popen(cmd, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=None, close_fds=True)
		proc_data = proc.communicate()[0]
		print proc_data
		data = 1
	
	except Exception, e:
		data = e
	return data

#for each moment concat moments for subject, generate temp matrix file
moments_list = ['mean','stdv','skew','kurt']

# svd with log(cosh) or svds with gibbs free energy
decomp = 'svds'

selected_moment = 0
people_dict = dict()
people_moment_dict = dict()

data_file_path = "filtered_obs"

moments_files = [x for x in os.listdir(data_file_path) if '.mat' in x]

for x in moments_files:
	subject = x[:2]

	if moments_list[selected_moment] + '.mat' in x:
		if subject not in people_dict:
			people_dict[subject] = [x]
		else:
			people_dict[subject].append(x)


for key, value in people_dict.iteritems():
	data_string = []
	out_concat_file = data_file_path + "/" + key + "/concat_" +moments_list[selected_moment] +".mat"

	if not isfile(out_concat_file):
		for x in value:
			temp_data = filter(None,open(data_file_path + "/" + x).read().split('\n'))
			#print len(temp_data)
			data_string += temp_data


	
	if not isdir(data_file_path + "/" + key):
		os.mkdir(data_file_path + "/" + key)

	if not isfile(out_concat_file):
		open(out_concat_file,'w+').write('\n'.join(data_string) + '\n')

	people_moment_dict[key] = out_concat_file

for key,value in people_moment_dict.iteritems():
	print 'Generating beamformer weights for subject %s' % (key)
	out_beaformer_weight_file = data_file_path + "/" + key + "/beamformer_weights_" +decomp +"_" +moments_list[selected_moment] +".mat"

	if not isfile(out_beaformer_weight_file):
		cmd1 = "./generate_covariance_matrix %s" % (value)
		cmd2 = "./forward_solution %s" % (out_beaformer_weight_file)
		#run generate_covariance_matrix with concated moment file
		r_cmd1 = execute_cmd(cmd1)

		if r_cmd1 == 1:
			#compute beamformer weights with save weights in folder for subject moment file and save weights in folder for subject
			execute_cmd(cmd2)
	
	sys.exit()