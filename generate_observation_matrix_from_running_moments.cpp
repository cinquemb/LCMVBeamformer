#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
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
			std::string temp_data_file = sub_file.asString();
			Json::Value moments_data = load_json(temp_data_file);
		}
	}

}
