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

int main(int argc, char *argv[]){
	//load in running moments files for each subject

	//load in binned for each data file into map<int,vector<double>>

	//iterate over subject and for each moment file for specified type of run, generate observation matricies and save files

}
