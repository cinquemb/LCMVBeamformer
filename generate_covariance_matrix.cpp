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
const int window_size = 2048;

//factor to convert biosemi values into uv
double biosemi_microvoltage_factor = 8192;

//num_colums
int num_colums = 128;

/*
# XXX: these should be optimized, as they can be a bottleneck.
def _logcosh(x, fun_args=None):
    alpha = fun_args.get('alpha', 1.0) # comment it out?
    x *= alpha
    gx = np.tanh(x, x) # apply the tanh inplace
    g_x = np.empty(x.shape[0])
    # XXX compute in chunks to avoid extra allocation
    for i, gx_i in enumerate(gx): # please don't vectorize.
        g_x[i] = (alpha * (1 - gx_i ** 2)).mean()
    return gx, g_x
*/

arma::mat logcosh(arma::mat& x){
    arma::mat gx = arma::tanh(x);
    return gx;
};

arma::mat sym_decorrelation(arma::mat& w){
    //W <- (W * W.T) ^{-1/2} * W
    arma::vec eigval;
    arma::mat eigvec;
    arma::mat w_dot = (w * w.t());
    std::cout << "post w_dot" <<std::endl;
    arma::eig_sym(eigval, eigvec, w_dot);
    std::cout << "post eig_sym" <<std::endl;
    arma::mat n_w = (((eigvec * arma::diagmat(1 / arma::sqrt(eigval))) * eigvec.t()) * w);
    return n_w;
}

arma::mat fast_ica_parallel(arma::mat& x, arma::mat& w_init, int max_iter, double tolerance){
    /*
    """Parallel FastICA.
Used internally by FastICA --main loop
"""
W = _sym_decorrelation(w_init)
del w_init
p_ = float(X.shape[1])
for ii in moves.xrange(max_iter):
    gwtx, g_wtx = g(fast_dot(W, X), fun_args)
    W1 = _sym_decorrelation(fast_dot(gwtx, X.T) / p_
        - g_wtx[:, np.newaxis] * W)
    del gwtx, g_wtx
    # builtin max, abs are faster than numpy counter parts.
    lim = max(abs(abs(np.diag(fast_dot(W1, W.T))) - 1))
    W = W1
    if lim < tol:
        break
else:
warnings.warn('FastICA did not converge. Consider increasing '
'tolerance or the maximum number of iterations.')
return W, ii + 1
    */
    arma::mat w = sym_decorrelation(w_init);
    std::cout << "post sym_decorrelation" <<std::endl;
    double p_ = w_init.n_cols;
    for(int i=0;i<max_iter;++i){
        arma::mat wx = w * x;
        std::cout << "post wx" <<std::endl;

        arma::mat gx = logcosh(wx);
        std::cout << gx <<std::endl;
        exit(0);

    }
}

arma::mat whiten_matrix_samples(arma::mat& sample_matrix){
    //column normalize matrix
    int n_components = std::min({sample_matrix.n_rows, sample_matrix.n_cols});
    arma::mat col_norm_sample_matrix = sample_matrix;
    for(int i=0;i<n_components;++i)
        col_norm_sample_matrix.row(i) = sample_matrix.row(i) - arma::mean(sample_matrix.row(i));

    arma::mat u;
    arma::vec simga_matrix;
    arma::mat v_matrix;
    arma::svds(u, simga_matrix, v_matrix, arma::sp_mat(col_norm_sample_matrix), n_components);

    arma::mat k(u.n_rows,u.n_cols);
    for(int i=0;i<n_components;++i)
        k.col(i) = sample_matrix.col(i)/simga_matrix(i);

    arma::mat whiten_matrix = (k * col_norm_sample_matrix) * std::sqrt(sample_matrix.n_cols);
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
    return sample_matrix;
}

int main(int argc, char *argv[]) {
	std::string file_name(argv[1]);
	arma::mat sample_matrix = matrix_from_file_samples(file_name);
    arma::mat sample_matrix_t = sample_matrix.t();
    double smm = arma::mean(arma::mean(sample_matrix));
    arma::mat sample_matrix_norm = sample_matrix - smm;
    arma::mat whiten_matrix = whiten_matrix_samples(sample_matrix_t);
    
    arma::mat w_init = arma::randn(sample_matrix.n_cols, sample_matrix.n_cols);
    arma::mat covmat = arma::cov(sample_matrix_norm);
    fast_ica_parallel(whiten_matrix, covmat, 200, 1e-04);

	covmat.save("test_cov.mat", arma::raw_ascii);
	return 0;
}