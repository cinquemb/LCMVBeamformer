import os
import json

data_dir_data_name = '/Users/cinquemb/Documents/Startups/GoBlue/Research-Code/NFBcorrelationAnalysis/eeg_raw_distrubution_plots/nfb_histogram_data_real/nfb_running_moments/'
temp_dir_data = os.listdir(data_dir_data_name)


moments_list = ['mean','stdv','skew','kurt']
outliars = ['10','14','20','26','32','48','blank','DAE']
exclude_list = ['23_RVIP_01_090814_1004_Run1_raw_bin_size_100', '23_10_EMG_01_090814_0954_Run1_raw_bin_size_100', '23_01_RS_01_090814_1019_Run1_raw_bin_size_100']
white_list = ['20_10_EMG_02_082214_1019_Run1_raw_bin_size_100','20_10_EMG_01_082214_1008_Run1_raw_bin_size_100']
explore_files_dict = {}

for x in temp_dir_data:
	if x != '.DS_Store':
		_file_ = data_dir_data_name+x
		if '_Run1_raw_bin_size_100.json' in x:
			subject = x.split('_')[0]
			if x in white_list:
				if subject not in explore_files_dict:
					explore_files_dict[subject]= [_file_]
				else:
					explore_files_dict[subject].append(_file_)
				#let thru
				print x
				continue
			elif x in exclude_list:
				#dont parsse
				#print x
				continue
			elif any(subject == y for y in outliars):
				#dont parsse
				continue
				#print x
				
			if subject not in explore_files_dict:
				explore_files_dict[subject]= [_file_]
			else:
				explore_files_dict[subject].append(_file_)

print json.dumps(explore_files_dict)