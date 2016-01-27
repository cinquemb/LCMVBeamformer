import os
import sys
import re
import csv
import codecs
import math
from os.path import isdir, isfile, join
from subprocess import Popen
import subprocess
import re
from PIL import Image
import datetime

"""
SCRIPT THAT ANALYZES IMAGES FOR REGION OF INTEREST AND COUNTS IF EVENT OCCURED IN SUBSET OF REGION IN ORDER OF OCURRANCE, NORMALIZES FOR RANDOMIZATION

ASSUMPTIONS:
	- one graph per photo
	- graph occurs at center of image
	- baseline marker to give start
	- hardcoded values for region of interest
"""

is_nested_dir = False

def find_start_occuerances(string, pattern):
	return [str(m.start()) for m in re.finditer('(?=' + pattern + ')', string)]

def round_format(nom,denom):
	try:
		if float(denom) > 0:
			percent = 100* round(float(nom)/float(denom), 4)
		else:
			percent = 0
	except Exception, e:
		percent = 0

	return percent

def execute_cmd(cmd):
	try:
		proc = Popen(cmd, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=None, close_fds=True)
		proc_data = proc.communicate()[0]
		data = 1
	
	except Exception, e:
		data = e
	return data

def find_x_val_border_left(width, height, pixel_of_red):
	for i in range(0,height):
		for j in range(0, width):
			temp_pix = im.getpixel((j,i))
			check_pix = [i for k, z in zip(temp_pix, pixel_of_red) if k == z]

			if len(check_pix) == len(pixel_of_red):
				return j

def find_y_val_center_init(width, height):
	width = width+1
	while True:
		prev_pix = dict()
		temp_y = None
		for i in range(0,height-1):
			tempix = ' '.join(map(str, im.getpixel((width,i))))
			if i > 0:
				if tempix not in prev_pix:
					temp_y = i
					return [width,temp_y]
			else:
				prev_pix[tempix] = 0

		width +=1

def count_green(width, height_min, height_max,calc_area_under_bin=False):
	pixel_of_green = (0, 255, 0, 255)
	green_check = ' '.join(map(str, pixel_of_green))
	count = 0
	width = width
	for i in range(height_min, height_max):
		tempix = ' '.join(map(str, im.getpixel((width,i))))
		
		if calc_area_under_bin:
			temp_crop_tup = (width,i,width,height_max)
			colors = im.crop(temp_crop_tup).getcolors()
			colors_in_crop = [' '.join(map(str, color[1])) for color in colors]

			if green_check not in colors_in_crop:
				break
			
			if tempix == green_check:
				count +=1
		else:
			if tempix == green_check:
				count +=1
				break

	return count

def get_substrings(string, size):
	substring = []
	for i in range(0,len(string)-size):
		substring.append(string[i:i+size])

	return substring

def get_inversion_file(photo_file_string, subject_num, run_no):
	len_of_interest = 16
	for i in range(0,len(header_files_list)):
		pic_string = header_screen_shot_time_dict[header_files_list[i]]['photo']
		if pic_string in photo_file_string:
			return header_files_list[i]

def is_inverted(photo_file_string, header_files_list, subject_num, run_no):
	numb_to_inversion = {
		'0': False,
		'1': True
	}
	
	inversion_file = get_inversion_file(photo_file_string, subject_num, run_no)
	#print photo_file_string, subject_num, run_no
	#print photo_file_string
	'''
	if inversion_file is None:
		print photo_file_string
		sys.exit()
	'''
	with open(inversion_file, 'r+') as f_ss:
		reader = csv.reader(f_ss, dialect='excel', delimiter=' ')
		for k, row in enumerate(reader):
			if 'inversionFlag_' in ' '.join(row):
				letter_val = row[-1]
				#print letter_val
				return [numb_to_inversion[letter_val], inversion_file]


dir_data_name = '/Users/cinquemb/Documents/Startups/GoBlue/Research-Code/png_headers_matched/'

dir_data = os.listdir(dir_data_name)
raw_data_file_list = []
events_file_list = []
stdin_files = []

if is_nested_dir:
	pass
else:
	for tem_file_ in dir_data:
		#print tem_file_
		temp_file_path = '%s%s' % (dir_data_name, tem_file_)
		stdin_files.append(temp_file_path)

query_of_interest = ['MCG','MGC']
raw_query_filt_data_file_list = []
header_files_list = []
for ss_file in stdin_files:
	#print ss_file
	if '.png' in ss_file and any(x in ss_file for x in query_of_interest):
		size = os.stat(ss_file).st_size
		#print size
		if size > 5*(10**3):
			raw_data_file_list.append(ss_file)
			raw_query_filt_data_file_list.append(ss_file)
	
	if 'HEADER' in ss_file:
		if 'MCG' in ss_file:
			header_files_list.append(ss_file)
		if 'MGC' in ss_file:
			header_files_list.append(ss_file)

header_screen_shot_time_dict = dict()
for x in header_files_list:
	h_d =  open(x, 'r+')
	t_h_d = filter(None,h_d.read().split('\n'))

	h_data = [i for i in  t_h_d if 'Screen_shot' in i][0].split(' ')[-1]
	h_r_data = [i for i in  t_h_d if 'Raw_Data' in i][0].split(' ')[-1].split('/')[-1]
	h_d.close()
	header_screen_shot_time_dict[x] = {'photo':h_data, 'data': h_r_data}
	

if not is_nested_dir:
	header_dirs = []
	d_dir_data = []
	for dirs_ in header_dirs:
		tmp_dirs = os.listdir(dirs_)
		out_dirs = []
		for path in tmp_dirs:
			temp_path_string = '%s%s/' % (dirs_, path)
			out_dirs.append(temp_path_string)

		d_dir_data.append(out_dirs)

	nested_header_dirs = [dirs for sub_node in d_dir_data for dirs in sub_node]
	stdin_hfiles = []

	for dir_ in nested_header_dirs:
		temp_dir_path = dir_
		if isdir(join(temp_dir_path)) and 'Expanded' not in temp_dir_path:
			temp_dir_data = os.listdir(temp_dir_path)
			for tem_file_ in temp_dir_data:
				#print temp_dir_path + '/' + tem_file_
				temp_dir_path_mod = '/'.join(temp_dir_path.split('/'))
				#print tem_file_
				temp_file_path = '%s%s' % (temp_dir_path_mod, tem_file_)
				stdin_hfiles.append(temp_file_path)

	for hh_file in stdin_hfiles:
		if 'HEADER' in hh_file:
			if 'MCG' in hh_file:
				header_files_list.append(hh_file)
			if 'MGC' in hh_file:
				header_files_list.append(hh_file)
 
filter_people = ['14','15','16','17','20','22','24','27','28','29','30','34','35','36','37', '38','40','41','42','43','45','46','47','48','50','51','54','55','57','58','59','61']

filter_people_filt_data_file_list = []

for interest_file in raw_query_filt_data_file_list:
	temp_out_file_string = (interest_file.split('.png')[0]).split('/')[-1] + '.png'
	if len(filter_people) > 0:
		if not any("ScreenShot__" + word == temp_out_file_string[0:14] for word in filter_people):
			continue
		else:
			filter_people_filt_data_file_list.append(interest_file)
	else:
		if int(temp_out_file_string[12:14]) < 14:
			continue
		else:
			filter_people_filt_data_file_list.append(interest_file)

i = 0

f = open('raw_files_list_with-%s-of-interest-pics.txt' % ('_'.join(query_of_interest)), 'w')
f.write('\n'.join(filter_people_filt_data_file_list))
f.close()

		
print 'Subject #, MCG1/2, Run #, %Direction Opposite EA, %Direction EA, Inversion'

out_raw_data_dict = dict()

for interest_file in filter_people_filt_data_file_list:

	temp_out_file_string = (interest_file.split('.png')[0]).split('/')[-1] + '.png'
	

	subject_num = temp_out_file_string[12:14]
	run_name = temp_out_file_string[18:22]
	run_no = temp_out_file_string[23:25]

	#print 'run_no', run_no

	if ('MCG' not in run_name) and ('MGC' not in run_name):
		if 'MCG1' in temp_out_file_string:
			start_offset = int(find_start_occuerances(temp_out_file_string, 'MCG1')[0]) + 5
			run_name = 'MCG1'
			run_no = temp_out_file_string[start_offset:start_offset+2]
			#print 'run_no mal', run_no
		elif 'MCG2' in temp_out_file_string:
			start_offset = int(find_start_occuerances(temp_out_file_string, 'MCG2')[0]) + 5
			run_name = 'MCG2'
			run_no = temp_out_file_string[start_offset:start_offset+2]
			#print 'run_no mal', run_no

	inversion, header_file = is_inverted(temp_out_file_string, header_files_list, subject_num, run_no)
	raw_file = header_screen_shot_time_dict[header_file]['data']


	#copycmd = 'cp %s /Users/cinquemb/Documents/Startups/GoBlue/Research-Code/interest_photos/' % (re.escape(interest_file))
	#execute_cmd(copycmd)
	'''
	replace with c++ code for area under curve using image magick? or let prasanta do it? 
	'''
	above = 0
	below = 0
	im = Image.open(interest_file)

	pixel_of_red = (255, 0, 0, 255)
	pixel_of_black = (0, 0, 0, 255)
	pixel_of_green = (0, 255, 0, 255)

	width = im.size[0]
	height = im.size[1]
	ofset_width = find_x_val_border_left(width, height, pixel_of_red)

	i = 0
	is_activated_bins = []
	activation_start_flag = False
	activation_start_val = 0
	split_region_val = find_y_val_center_init(ofset_width, height)

	for j in range(split_region_val[0], width):
		temp_crop_tup = (j,0,width,height)
		green_exist_check = ' '.join(map(str, pixel_of_green))
		colors = im.crop(temp_crop_tup).getcolors()
		colors_in_crop = [' '.join(map(str, color[1])) for color in colors]

		if green_exist_check not in colors_in_crop:
			break

		green_count_below = count_green(j, height/2, height)
		green_count_above = count_green(j, 0, height/2)
		if green_count_above > 0:
			if not activation_start_flag:
				activation_start_flag = True
				activation_start_val = i
			above += 1

		if green_count_below > 0:
			if activation_start_flag:
				is_activated_bins.append([activation_start_val, i])
				activation_start_flag = False
			below += 1

		i+=1

	total = above + below
	per_above = round_format(above, total)
	per_below = round_format(below, total)

	if inversion:
		print 'Deactivation'
		temp_swap_b = per_above
		temp_swap_a = per_below
		per_below = temp_swap_b
		per_above = temp_swap_a
	else:
		print 'Activation'

	#print 'above: ',above, 'below: ',below, 'percent above: ',per_above,'percent below: ',per_below ,re.escape(interest_file)
	#print '%s, %s, %s, %s, %s, %s' % (subject_num, run_name, run_no, per_above, per_below, str(inversion))
	out_raw_data_dict[raw_file] = {
		'subject_num': subject_num,
		'run_name': run_name,
		'run_no': run_no,
		'per_above': per_above,
		'per_below': per_below,
		'inversion': inversion,
		'inversion_canoical': 'Deactivation' if inversion else 'Activation',
		'header_file': header_file,
		'photo_file': interest_file,
		'is_activated_bin_sans_norm': is_activated_bins
	}
	#print out_raw_data_dict[raw_file]
	#print '\n\n\n\n'
	print interest_file, 'Done'

open('mined_activation_markers.json','w+').write(json.dumps(out_raw_data_dict))