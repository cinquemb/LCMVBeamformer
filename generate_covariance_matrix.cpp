#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <chrono>
#include <armadillo>

#define IS_TIMER_ON

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
    arma::mat u, v_matrix; 
    arma::vec simga_matrix;
    arma::svds(u, simga_matrix, v_matrix, arma::sp_mat(w), n_components);
    return (u * v_matrix.t());
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
    int n_components = std::min({sample_matrix.n_rows, sample_matrix.n_cols});
    arma::mat col_norm_sample_matrix = sample_matrix;

    arma::mat u, v_matrix;
    arma::vec simga_matrix;
    arma::svds(u, simga_matrix, v_matrix, arma::sp_mat(col_norm_sample_matrix), n_components);

    arma::mat whiten_matrix = u * v_matrix.t();
    return whiten_matrix;
}

int main(int argc, char *argv[]) {
	std::string file_name(argv[1]);
	arma::mat sample_matrix;

    sample_matrix.load(file_name, arma::raw_ascii);

    #ifdef IS_TIMER_ON
    auto normalization_start = std::chrono::high_resolution_clock::now();
    #endif
    double smmin = arma::min(arma::min(sample_matrix));
    double smmax = arma::max(arma::max(sample_matrix));
    double sm_range = smmax - smmin;

    arma::mat sample_matrix_norm = (sample_matrix - smmin)/sm_range;
    #ifdef IS_TIMER_ON
    auto normalization_end = std::chrono::high_resolution_clock::now();
    auto normalization_difftime = normalization_end - normalization_start;
    std::cout << "normalization difftime (μs): " << std::chrono::duration_cast<std::chrono::microseconds>(normalization_difftime).count() << std::endl;
    #endif

    arma::mat sample_matrix_norm_t = sample_matrix_norm.t();
    #ifdef IS_TIMER_ON
    auto whiten_start = std::chrono::high_resolution_clock::now();
    #endif
    arma::mat whiten_matrix = whiten_matrix_samples(sample_matrix_norm_t);
    #ifdef IS_TIMER_ON
    auto whiten_end = std::chrono::high_resolution_clock::now();
    auto whiten_difftime = whiten_end - whiten_start;
    std::cout << "whitten difftime (μs): " << std::chrono::duration_cast<std::chrono::microseconds>(whiten_difftime).count() << std::endl;
    #endif
    
    arma::mat w_init = arma::randn(sample_matrix.n_cols, sample_matrix.n_cols);
    #ifdef IS_TIMER_ON
    auto covmat_start = std::chrono::high_resolution_clock::now();
    #endif
    arma::mat covmat = arma::cov(sample_matrix_norm);
    #ifdef IS_TIMER_ON
    auto covmat_end = std::chrono::high_resolution_clock::now();
    auto covmat_difftime = covmat_end - covmat_start;
    std::cout << "covmat difftime (μs): " << std::chrono::duration_cast<std::chrono::microseconds>(covmat_difftime).count() << std::endl;
    #endif
    
    #ifdef IS_TIMER_ON
    auto start = std::chrono::high_resolution_clock::now();
    #endif
    arma::mat ica_matrix = fast_ica_parallel(whiten_matrix, w_init, 200, 5e-2);
    #ifdef IS_TIMER_ON
    auto end = std::chrono::high_resolution_clock::now();
    auto difftime = end - start;
    std::cout << "fast_ica_parallel difftime (μs): " << std::chrono::duration_cast<std::chrono::microseconds>(difftime).count() << std::endl;
    #endif

	covmat.save("test_cov.mat", arma::raw_ascii);
    ica_matrix.save("ica_matrix.mat", arma::raw_ascii);
	return 0;
}