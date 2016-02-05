#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <iterator>
#include <map>

#include <armadillo>

#include <json/json.h>
#include <json/reader.h>

#include <boost/algorithm/string.hpp>

const int offset_node = 0;

//window size
const int window_size = 2048;

//factor to convert biosemi values into uv
double biosemi_microvoltage_factor = 8192;

//num_colums
int num_colums = 128;

Json::Value load_json(std::string file){
    Json::Value root;
    {
        std::ifstream json_data(file);
        try{
            json_data >> root;
        }catch(...){
             Json::Value root;
        }
    }
    return root;
}

int main(int argc, char *argv[]){
	//load in running moments files for each subject
	Json::Value running_moments_binned_series = load_json("running_moments_ts.json");	

	//load in binned for each data file into map<int,vector<double>>
	Json::Value activation_markers = load_json("mined_activation_markers.json");
	
	//iterate over subject and for each moment file for specified type of run, generate observation matricies and save files
	for(auto sub_id : running_moments_binned_series){
		for(auto sub_file : sub_id){
			std::string temp_moment_data_file = sub_file.asString();

			/*

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

			*/
			std::vector<std::string> temp_data_file_strings;
			boost::split(temp_data_file_strings, temp_data_file, boost::is_any_of(","), boost::token_compress_on);
			std::string temp_raw_data_file = temp_data_file_strings[temp_data_file_strings.size()-1];
			Json::Value temp_activation_marker_data = activation_markers[temp_raw_data_file];

			std::string temp_inversion_canoical = temp_activation_marker_data['inversion_canoical'].asString();

			if(temp_inversion_canoical == 'Activation'){
				//generate inverse marked bins and save as marked bins in vector of vectors (inner vector is [start index, end index])
				//expand marked bins into time bins

			}else{
				//expand marked bins into time bins
			}

			Json::Value moments_data = load_json(temp_moment_data_file);
		}
	}

}
