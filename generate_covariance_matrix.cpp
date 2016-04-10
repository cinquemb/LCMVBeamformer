#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <map>
#include <chrono>

#include <armadillo>

#include <boost/algorithm/string.hpp>

const int offset_node = 0;

//window size
const int window_size = 2048;

//factor to convert biosemi values into uv
double biosemi_microvoltage_factor = 8192;

//num_colums
int num_colums = 128;

bool test_mode = true;

bool matrix_from_moment = true;

std::vector<arma::mat> logcosh(arma::mat& x){
    arma::mat gx = arma::tanh(x);
    arma::mat g_x(gx.n_rows, 1);
    for(int i=0; i<gx.n_rows; ++i){
        arma::rowvec dtanh_i = 1 -arma::pow(gx.row(i),2);
        g_x.row(i) = arma::mean(dtanh_i);
    }
    std::vector<arma::mat> gx_dgx{gx,g_x};
    return gx_dgx;
};

arma::mat sym_decorrelation(arma::mat& w){
    //W <- (W * W.T) ^{-1/2} * W or W = U * V.t()
    int n_components = std::min({w.n_rows, w.n_cols});
    arma::mat u;
    arma::vec simga_matrix;
    arma::mat v_matrix;

    if(!test_mode)
        arma::svd(u, simga_matrix, v_matrix, w);
    else
        arma::svds(u, simga_matrix, v_matrix, arma::sp_mat(w), n_components);

    return (u * v_matrix.t());
}

arma::mat sym_decorrelation_complex(arma::mat& w){
    //W <- (W * W.T) ^{-1/2} * W

    arma::cx_mat ww = arma::cx_mat(w,arma::mat(w.n_rows, w.n_cols, arma::fill::zeros));
    arma::vec eigval;
    arma::cx_mat eigvec;
    arma::cx_mat w_dot = (ww * ww.t());
    arma::eig_sym(eigval, eigvec, w_dot);
    arma::cx_mat n_w = (((eigvec * arma::diagmat(1 / arma::sqrt(eigval))) * eigvec.t()) * ww);
    arma::cx_mat n_w_h_adj = arma::conj(n_w.t());
    arma::cx_mat r_n_w = n_w * n_w_h_adj;
    return arma::real(r_n_w);
}

arma::mat fast_ica_parallel(arma::mat& x, arma::mat& w_init, int max_iter, double tolerance){
    arma::mat w = sym_decorrelation(w_init);
    double p_ = w_init.n_cols;
    for(int i=0;i<max_iter;++i){
        arma::mat wx = w * x;
        std::vector<arma::mat> gx_dgx = logcosh(wx);
        arma::mat gwtx_xt = (gx_dgx[0] * x.t()) / p_;
        
        arma::mat g_wtx_w = arma::repmat(gx_dgx[1], 1, gx_dgx[1].n_rows) * w;
        arma::mat t_w1 = gwtx_xt - g_wtx_w;

        arma::mat w1 = sym_decorrelation(t_w1);
        //std::cout << t_w1 << std::endl;

        arma::mat w1_wt = w1 * w.t();
        double lim = arma::max(arma::abs((arma::abs((w1_wt.diag())) -1)));
        std::cout << "lim: " << lim << " tolerance: "  << tolerance << std::endl;
        if(lim < tolerance)
            break;

        if(i == max_iter-1)
            std::cout << "max iteration excceeded" << std::endl;
        w = w1;   
    }
    return w;
}

arma::mat whiten_matrix_samples(arma::mat& sample_matrix){
    //column normalize matrix
    int n_components = std::min({sample_matrix.n_rows, sample_matrix.n_cols});
    arma::mat col_norm_sample_matrix = sample_matrix;

    if(!test_mode){
        for(int i=0;i<n_components;++i)
            col_norm_sample_matrix.row(i) = sample_matrix.row(i) - arma::mean(sample_matrix.row(i));
    }else{
        for(int i=0;i<sample_matrix.n_cols;++i){
            double col_min = arma::min(sample_matrix.col(i));
            double col_max = arma::max(sample_matrix.col(i));
            double col_diff = col_max - col_min;
            col_norm_sample_matrix.col(i) = (sample_matrix.col(i) - col_min)/col_diff;
        }
            
    }
    

    arma::mat u;
    arma::vec simga_matrix;
    arma::mat v_matrix;

    if(!test_mode)
        arma::svd(u, simga_matrix, v_matrix, col_norm_sample_matrix);
    else
        arma::svds(u, simga_matrix, v_matrix, arma::sp_mat(col_norm_sample_matrix), n_components);
    /*
    arma::mat k(u.n_rows,u.n_cols);
    for(int i=0;i<n_components;++i)
        k.col(i) = u.col(i)/simga_matrix(i);

    arma::mat whiten_matrix = (k * col_norm_sample_matrix) * std::sqrt(sample_matrix.n_cols);
    */
    arma::mat whiten_matrix = u * v_matrix.t();
    return whiten_matrix;
}

arma::mat matrix_from_file_samples(std::string& file_name){
	std::vector<std::string> data_vector_string;
	std::vector<double> data_vector_out;
	std::string line;

	int line_count = 0;

	std::ifstream in(file_name.c_str());
	if (!in.is_open()) exit(0);

	arma::mat sample_matrix(window_size, num_colums);

	while (std::getline(in,line)){
		if(line_count >= window_size)
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
    return sample_matrix;
}

int main(int argc, char *argv[]) {
    if(test_mode)
        std::cout << "test mode" << std::endl;
	std::string file_name(argv[1]);
	arma::mat sample_matrix;

    if(matrix_from_moment)
        sample_matrix.load(file_name, arma::raw_ascii);
    else
        sample_matrix = matrix_from_file_samples(file_name);

    auto init_start = std::chrono::high_resolution_clock::now();
    arma::mat sample_matrix_t = sample_matrix.t();
    double smmin = arma::min(arma::min(sample_matrix));
    double smmax = arma::max(arma::max(sample_matrix));
    double sm_range = smmax - smmin;

    arma::mat sample_matrix_norm = (sample_matrix - smmin)/sm_range;


    arma::mat sample_matrix_norm_t = sample_matrix_norm.t();
    arma::mat whiten_matrix = whiten_matrix_samples(sample_matrix_norm_t);
    
    arma::mat w_init = arma::randn(sample_matrix.n_cols, sample_matrix.n_cols);
    arma::mat covmat = arma::cov(sample_matrix_norm);
    auto init_end = std::chrono::high_resolution_clock::now();
    auto init_difftime = init_end - init_start;
    std::cout << "init preprocess difftime (μs): " << std::chrono::duration_cast<std::chrono::microseconds>(init_difftime).count() << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    arma::mat ica_matrix = fast_ica_parallel(whiten_matrix, w_init, 200, 5e-2);
    auto end = std::chrono::high_resolution_clock::now();
    auto difftime = end - start;
    std::cout << "fast_ica_parallel difftime (μs): " << std::chrono::duration_cast<std::chrono::microseconds>(difftime).count() << std::endl;

	covmat.save("test_cov.mat", arma::raw_ascii);
    ica_matrix.save("ica_matrix.mat", arma::raw_ascii);
	return 0;
}