import os
import sys
from os.path import isdir, isfile, join
from subprocess import Popen
import subprocess
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats.stats import pearsonr
import numpy as np

#find . -maxdepth 3 -type f -name '*svd_*mean*.png' -print0 | xargs -0n 100000 open
#find . -maxdepth 3 -type f -name '*svd_*stdv*.png' -print0 | xargs -0n 100000 open
#find . -maxdepth 3 -type f -name '*svd_*skew*.png' -print0 | xargs -0n 100000 open
#find . -maxdepth 3 -type f -name '*svd_*kurt*.png' -print0 | xargs -0n 100000 open

#find . -maxdepth 3 -type f -name '*svd_whitten_diff_5*_mean*.png' -print0 | xargs -0n 100000 open
#find . -maxdepth 3 -type f -name '*svd_whitten_diff_5*_stdv*.png' -print0 | xargs -0n 100000 open
#find . -maxdepth 3 -type f -name '*svd_whitten_diff_5*_skew*.png' -print0 | xargs -0n 100000 open
#find . -maxdepth 3 -type f -name '*svd_whitten_diff_5*_kurt*.png' -print0 | xargs -0n 100000 open

#for each moment concat moments for subject, generate temp matrix file
moments_list = ['mean','stdv','skew','kurt']

# svd with log(cosh) or svds with gibbs free energy
decomp = 'svd_whitten_diff_5'

def get_stdv(values, mean):
	dif_sqaured = []
	for val in values:
		dfsq = (val-mean)**2
		dif_sqaured.append(dfsq)

	stdV = math.sqrt(sum(dif_sqaured)/float(len(dif_sqaured)-1))

	return stdV

def get_correlation(data_x, data_y):
	mu = 2
	len_data_x = len(data_x)
	len_data_y = len(data_y)
	if len_data_x != len_data_y:
		if len_data_y > len_data_x:
			data_y = data_y[:len_data_x]
		else:
			data_x = data_x[:len_data_y]

	mean_x = sum(data_x)/float(len(data_x))
	mean_y = sum(data_y)/float(len(data_y))

	data_x_stdv = get_stdv(data_x, mean_x)
	data_y_stdv = get_stdv(data_y, mean_y)
	
	outliars_index = []
	data_y_normalized_outliars_index = []
	for i,x in enumerate(data_x):
		if (abs(x-mean_x) > (mu*data_x_stdv)) or (abs(data_y[i]-mean_y) > (mu*data_y_stdv)):
			outliars_index.append(i)

	#print outliars_index
	x1 = [x for i,x in enumerate(data_x) if i not in outliars_index]
	y1 = [y for i,y in enumerate(data_y) if i not in outliars_index]
	#print len(x1), len(y1)

	correlation = pearsonr(data_x[:-3],data_y[:-3])[0]
	#print correlation,data_x_stdv, data_y_stdv
	return correlation

def execute_cmd(cmd):
	try:
		proc = Popen(cmd, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=None, close_fds=True)
		proc_data = proc.communicate()[0]
		print proc_data
		data = 1
	
	except Exception, e:
		data = e
	return data


people_dict = dict()
people_moment_dict = dict()

data_file_path = "filtered_obs"
exp_bwf = "experimental_beamformer_weights.txt"


moments_files = [x for x in os.listdir(data_file_path) if '.mat' in x]

for m_key in range(0, len(moments_list)):
	selected_moment = m_key
	print moments_list[selected_moment]
	for x in moments_files:
		subject = x[:2]
		if subject != '27':
			pass#continue
		else:
			pass

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
				temp_data = filter(None,open(data_file_path + "/" + x,'r+').read().split('\n'))
				#print len(temp_data)
				data_string += temp_data


		
		if not isdir(data_file_path + "/" + key):
			os.mkdir(data_file_path + "/" + key)

		if not isfile(out_concat_file):
			open(out_concat_file,'w+').write('\n'.join(data_string) + '\n')

		people_moment_dict[key] = out_concat_file

	out_beaformer_weight_file_map = dict()
	for key,value in people_moment_dict.iteritems():
		print 'Generating beamformer weights for subject %s' % (key)
		out_beaformer_weight_file = data_file_path + "/" + key + "/beamformer_weights_" +decomp +"_" +moments_list[selected_moment] +".mat"
		out_beaformer_weight_file_map[key] = out_beaformer_weight_file

		if not isfile(out_beaformer_weight_file):
			cmd1 = "./generate_covariance_matrix %s" % (value)
			cmd2 = "./forward_solution %s" % (out_beaformer_weight_file)
			#run generate_covariance_matrix with concated moment file
			r_cmd1 = execute_cmd(cmd1)

			if r_cmd1 == 1:
				#compute beamformer weights with save weights in folder for subject moment file and save weights in folder for subject
				execute_cmd(cmd2)
			else:
				print 'Failed: %s' % (out_beaformer_weight_file)
				sys.exit()


	exp_bwf_data = filter(None,open(exp_bwf).read().split('\n'))
	temp_sum = {0:[],1:[],2:[]}
	file_count = 0
	for key, value in out_beaformer_weight_file_map.iteritems():
		sub_exp_bwf_data = filter(None,open(value).read().split('\n'))
		print '\n\n'
		print 'Generating plot for', key
		for dim in range(0,3):
			temp_dim_exp = map(float, exp_bwf_data[dim].split('\t'))
			temp_sub_exp_bwf_data = filter(None,sub_exp_bwf_data[dim].split(' '))
			temp_dim_sub_exp = map(float,temp_sub_exp_bwf_data)
			chan_correlation = get_correlation(list(abs(np.array((temp_dim_exp)))), list(abs(np.array((temp_dim_sub_exp)))))

			tmp_dim = temp_sum[dim]
			temp_sum[dim] = list(np.sum([tmp_dim,temp_dim_sub_exp], axis=0))
			print 'correlation for channel %s: %s' % (dim, chan_correlation)
			out_beaformer_weight_file_png = data_file_path + "/" + key + "/" + key +"_beamformer_weights_" +decomp +"_" +moments_list[selected_moment] +"_dim_"+ str(dim) + ".png"

			if not isfile(out_beaformer_weight_file_png):
				fig = plt.figure(num=None, figsize=(25, 15), dpi=20, facecolor='w', edgecolor='k')
				plt.plot([x for x in range(len(temp_dim_exp))],temp_dim_exp, hold=True, color='red')
				plt.plot([x for x in range(len(temp_dim_exp))], temp_dim_sub_exp, color='blue')
				plt.ylim(-.05,.05)
				plt.xlim(0,128)
				plt.xlabel("lead index")
				plt.ylabel("weight values")
				fig.savefig(out_beaformer_weight_file_png)
				plt.close()

		file_count += 1
		print 'Plot created for', key
		print '\n\n'

	for key,value in temp_sum.iteritems():
		temp_dim_exp = map(float, exp_bwf_data[key].split('\t'))
		tmp_data = list(value)
		for x in range(0,len(tmp_data)):
			tmp_data[x] = tmp_data[x]/float(file_count)
		chan_correlation = get_correlation(list(abs(np.array((temp_dim_exp)))), list(abs(np.array((tmp_data)))))
		print 'avg correlation for channel %s: %s' % (key, chan_correlation)
