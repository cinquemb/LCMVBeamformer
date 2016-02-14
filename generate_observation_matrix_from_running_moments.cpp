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
				60 second feedbak window, 4 px per window so 240 max
				skip first 299 indecies for each channel
				
			'subject_num': subject_num,
			'run_name': run_name,
			'run_no': run_no,
			'per_above': per_above,
			'per_below': per_below,
			'inversion': inversion,
			'inversion_canoical': 'Deactivation' if inversion else 'Activation',
			'header_file': header_file,
			'photo_file': interest_file,
			'is_activated_bin_sans_norm': is_activated_bins,
			'is_deactivated_bin_sans_norm': is_deactivated_bin

			*/
			std::vector<std::string> temp_data_file_strings;

			boost::split(temp_data_file_strings, temp_moment_data_file, boost::is_any_of("/"), boost::token_compress_on);
			std::string temp_raw_data_file = temp_data_file_strings[temp_data_file_strings.size()-1];

			boost::split(temp_data_file_strings, temp_raw_data_file, boost::is_any_of("_"), boost::token_compress_on);
			temp_raw_data_file = "";
			for(int i =0; i< temp_data_file_strings.size(); ++i){
				if(temp_data_file_strings[i] == "bin"){
					temp_raw_data_file = temp_raw_data_file.substr(0, temp_raw_data_file.size()-1);
					break;
				}else
					temp_raw_data_file += temp_data_file_strings[i];

				temp_raw_data_file += "_";
			}

			temp_raw_data_file += ".txt";
			

			Json::Value temp_activation_marker_data = activation_markers.get(temp_raw_data_file,false);
			if(temp_activation_marker_data == false)
				continue;
			else
				std::cout << temp_raw_data_file << std::endl;


			bool temp_inversion_canoical_flag = temp_activation_marker_data["inversion"].asBool();

			Json::Value deactivated_bins;

			if(temp_inversion_canoical_flag)
				deactivated_bins = temp_activation_marker_data["is_activated_bin_sans_norm"];
			else
				deactivated_bins = temp_activation_marker_data["is_deactivated_bin_sans_norm"];

			/*
				for each window, 
					- add one to the endpoint
					- divide init by 4 to get offset index, normalize offset to 100ms bins
					- subtract endpoint from begining
					
					- divide diff by 4 to get total seconds, divide by .10 to get 100ms bins, add to offset to get end index
					- push back map of start and end 
			*/

			std::vector<std::map<std::string, int>> roi_time_series_boundaries;
			for(auto deac_pix_index : deactivated_bins){
				int mod_start = deac_pix_index[0].asInt()/4/.10;
				int mod_endpoint = (deac_pix_index[1].asInt() + 1)/4/.10;
				int mod_diff = mod_endpoint - mod_start;
				int real_endpint = mod_start + mod_diff;

				std::map<std::string, int> tmp_boundaries;
				tmp_boundaries["start_iter"] = mod_start;
				tmp_boundaries["end_iter"] = real_endpint;
				//std::cout << deac_pix_index[0] << "-" << deac_pix_index[1] << "|" << mod_start << "-" << real_endpint << std::endl;
				roi_time_series_boundaries.push_back(tmp_boundaries);
			}

			Json::Value moments_data = load_json(temp_moment_data_file);
			//for each channel
			for(Json::ValueIterator channel_id = moments_data.begin(); channel_id != moments_data.end(); channel_id++){
				int tmp_chan = atoi(channel_id.key())-2;
				std::map<std::string, arma::mat> moments_deactivation_samples;
				//for each moment for a given channel
				for(auto chan_moment : channel_id){
					//for each column value in each momenet for a given channel
					for(int i=0; i < chan_moment.size(); i++){
						//figure out how to properly insert values into sample matrix
						//moments_deactivation_samples[tmp_chan][i] = chan_moment[i];
					}
				}
			}
			/*
				iterate over each channel
					- for each moment in file, push back samples from start and end index into sample matrix
			*/
		}
	}

}
