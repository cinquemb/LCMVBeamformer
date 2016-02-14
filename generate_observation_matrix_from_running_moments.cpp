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

int offest_row_index = 300;

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

	std::vector<std::string> moments_list = {"mean","stdv","skew","kurt"};

	//load in binned for each data file into map<int,vector<double>>
	Json::Value activation_markers = load_json("mined_activation_markers.json");

	std::map<std::string, arma::mat> moments_deactivation_samples;
	
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
			std::string temp_raw_data_file_name = temp_raw_data_file;
			temp_raw_data_file += ".txt";
			

			Json::Value temp_activation_marker_data = activation_markers.get(temp_raw_data_file,false);
			if(temp_activation_marker_data == false)
				continue;
			else
				std::cout << "Starting mining: " << temp_raw_data_file << std::endl;


			bool temp_inversion_canoical_flag = temp_activation_marker_data["inversion"].asBool();

			Json::Value deactivated_bins;

			if(temp_inversion_canoical_flag)
				deactivated_bins = temp_activation_marker_data["is_activated_bin_sans_norm"];
			else
				deactivated_bins = temp_activation_marker_data["is_deactivated_bin_sans_norm"];

			std::vector<std::map<std::string, int>> roi_time_series_boundaries;
			for(auto deac_pix_index : deactivated_bins){
				int mod_start = deac_pix_index[0].asInt()/4/.10;
				int mod_endpoint = (deac_pix_index[1].asInt() + 1)/4/.10;
				int mod_diff = mod_endpoint - mod_start;
				int real_endpint = mod_start + mod_diff;

				std::map<std::string, int> tmp_boundaries;
				tmp_boundaries["start_iter"] = mod_start;
				tmp_boundaries["end_iter"] = real_endpint;
				roi_time_series_boundaries.push_back(tmp_boundaries);
			}

			Json::Value moments_data = load_json(temp_moment_data_file);
			std::map<std::string, arma::mat> tmp_moments_deactivation_samples;
			for(Json::ValueIterator channel_id = moments_data.begin(); channel_id != moments_data.end(); channel_id++){
				int tmp_chan = atoi(channel_id.key().asString().c_str())-2;
				for(auto tmp_chan_moment : moments_list){
					auto tmp_chan_moment_data = channel_id->get(tmp_chan_moment,false);
					if(tmp_moments_deactivation_samples.count(tmp_chan_moment) == 0){
						tmp_moments_deactivation_samples[tmp_chan_moment] = arma::mat(tmp_chan_moment_data.size(), num_colums);
					}

					for(int i=0; i < tmp_chan_moment_data.size(); i++){
						tmp_moments_deactivation_samples[tmp_chan_moment](i,tmp_chan) = tmp_chan_moment_data[i].asFloat();
					}
				}
			}

			std::cout << "Moment matricies loaded, starting slicing/filtering" << std::endl;
			std::map<std::string, arma::mat> tmp_moments_deactivation_samples_filtered;			
			for(auto moment_matrix : tmp_moments_deactivation_samples){
				arma::mat tmp_moment_deactivation_samples_filtered;
				std::cout << moment_matrix.second.n_rows << " by " << moment_matrix.second.n_cols << std::endl;
				int total_samples = 0;
				for(int i=0;i<roi_time_series_boundaries.size();++i){
					total_samples += roi_time_series_boundaries[i]["end_iter"] - roi_time_series_boundaries[i]["start_iter"];
					std::cout << "start: " << offest_row_index+roi_time_series_boundaries[i]["start_iter"] << "end: "<< offest_row_index+roi_time_series_boundaries[i]["end_iter"] << std::endl;
					arma::mat sub_tmp_moment_deactivation_samples_filtered = moment_matrix.second.rows(offest_row_index + roi_time_series_boundaries[i]["start_iter"], offest_row_index + roi_time_series_boundaries[i]["end_iter"]-1);
					arma::mat tmp_mdsf = arma::join_cols(tmp_moment_deactivation_samples_filtered, sub_tmp_moment_deactivation_samples_filtered);
					tmp_moment_deactivation_samples_filtered = tmp_mdsf;
					
				}
				tmp_moments_deactivation_samples_filtered[moment_matrix.first] = tmp_moment_deactivation_samples_filtered;
				std::cout << tmp_moment_deactivation_samples_filtered.n_rows << std::endl;
				std::cout << "total_samples: "<< total_samples << std::endl;
			}
			tmp_moments_deactivation_samples.clear();

			std::cout << "Done slicing/filtering, saving data" << std::endl;
			for(auto moment_matrix : tmp_moments_deactivation_samples_filtered){
				moment_matrix.second.save("filtered_obs/"+ temp_raw_data_file_name + "_"+ moment_matrix.first + ".mat", arma::raw_ascii);
				moment_matrix.second.clear();
			}
			std::cout << "Data saved\n\n\n" << std::endl;
		}
	}
}
