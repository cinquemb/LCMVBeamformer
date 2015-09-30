#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <map>

#include <armadillo>

#include <boost/algorithm/string.hpp>

const int offset_node = 0;

//window size
const int window_size = 2048*100;

//factor to convert biosemi values into uv
double biosemi_microvoltage_factor = 8192;

//num_colums
int num_colums = 128;

arma::mat avg_covariance_matrix_from_file_samples(std::string& file_name){
	std::vector<std::string> data_vector_string;
	std::vector<double> data_vector_out;
	std::string line;

	int line_count = 0;

	std::ifstream in(file_name.c_str());
	if (!in.is_open()) exit(0);

	arma::mat sample_matrix(window_size, num_colums);

	while (std::getline(in,line)){
		if(line_count >= 2048)
			break;

    	data_vector_string.clear();
    	data_vector_out.clear();

    	boost::split(data_vector_string,line,boost::is_any_of(" "));
    	for(int i =offset_node; i<data_vector_string.size();i++){
    		if (data_vector_string[i].size() > 0){
                double node_val;
                std::stringstream(data_vector_string[i]) >> node_val;
                data_vector_out.push_back(node_val/biosemi_microvoltage_factor);
            }      
    	}

    	std::vector<double> time_sample_data(data_vector_out.begin(), data_vector_out.begin()+num_colums);

    	sample_matrix.row(line_count) = arma::rowvec(time_sample_data);
    	line_count++;
    }
    in.close();
    arma::mat covariance_sample_matrix = arma::cov(sample_matrix);
    return covariance_sample_matrix;
}

int main(int argc, char *argv[]) {
	std::string file_name(argv[1]);
	arma::mat covmat = avg_covariance_matrix_from_file_samples(file_name);
	covmat.save("test_cov.mat", arma::raw_ascii);
	return 0;
}