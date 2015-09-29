#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <map>

#include <armadillo>

/*
REF:: http://cda.psych.uiuc.edu/matlab_class/html/5Toolbox/eeg_sph.html

INPUTS (Required):
       Rq   : dipole location(in meters)                                 P x 3 (variable)
       Re   : EEG sensors(in meters) on the scalp                        M x 3 (use biosemi map)

       R    : radii(in meters) of sphere from 
              INNERMOST to OUTERMOST                                     NL x 1
       sigma: conductivity from INNERMOST to OUTERMOST                   NL x 1

 INPUTS (Optional):
       nmax : # of terms used in Truncated Legendre Series               scalar
              If not specified, a default value based on outermost
              dipole magnitude is computed. (Note: This parameter
              is ignored when Berg Approximation is commanded)
       method : Method used for computing forward potential    
              1=Legendre Series Approx; 2=Berg Parameter Approx
              (Note: Default and all other values invoke Legendre 
              Series Approx. Exception is single-shell case where
              closed form solution is always used)                       scalar 
       mu_berg_in: User specified initial value for Berg eccentricity
              factors (Required if Berg Method is commanded)             3 x 1
	   lam_berg_in: User specified initial value for Berg magnitude
*/

const long double PI = acos(-1.0L);

double mill_to_meter = 1000.00;
double cent_to_meter = 100.00;
/* forward solution Eccentricy factors for 3 shell (mu_berg_in) */
std::vector<double> mu_berg = {.4191, .7479, .9791};

/*	magnitude factors (lam_berg_in) */
std::vector<double> lam_berg = {.4127, .2669, .0578};

/*  conductivities of regions (sigma) */
std::vector<double> sigmas = {.33, .0042, .33};

/*	radii fromm innermost to outer (R) */
std::vector<double> radii = {(71.00), (79.00), (85.00)};



/* eeg sensors in x,y,z using map (RE) */
std::vector<std::vector<double> > sensor_map =  {{-0.0,0.0,85.0}, {16.9462744255,1.03766003663e-15,83.2935998928}, {33.2121459216,2.03365740978e-15,78.2429125435}, {48.1445301386,2.94800223654e-15,70.0507260329}, {56.4895820686,-23.3987510256,59.045961489}, {43.2352543177,-43.2352543177,59.045961489}, {50.691263895,-50.691263895,45.6704667095}, {64.1990066054,-46.6433085545,30.4612757114}, {67.8233541008,-49.2765511462,14.0290464982}, {68.7245538619,-49.9313110988,-2.96645721971}, {66.8664221858,-48.5812994135,-19.8428559278}, {62.323564157,-45.2807198569,-35.922552248}, {73.2657437636,-23.8054832079,-35.922552248}, {78.6061936687,-25.5407005689,-19.8428559278}, {80.7905584608,-26.2504437136,-2.96645721971}, {79.731134603,-25.9062160355,14.0290464982}, {75.470458589,-24.5218384793,30.4612757114}, {66.2313280479,-27.4339143314,45.6704667095}, {61.1438830288,3.74398303193e-15,59.045961489}, {71.6882728941,4.38964069681e-15,45.6704667095}, {79.3543362523,4.85905169449e-15,30.4612757114}, {83.8342761307,5.13336889611e-15,14.0290464982}, {84.9482202966,5.20157830398e-15,-2.96645721971}, {82.6514432338,5.06094127006e-15,-19.8428559278}, {77.0361618981,4.71710445436e-15,-35.922552248}, {73.2657437636,23.8054832079,-35.922552248}, {78.6061936687,25.5407005689,-19.8428559278}, {80.7905584608,26.2504437136,-2.96645721971}, {79.731134603,25.9062160355,14.0290464982}, {75.470458589,24.5218384793,30.4612757114}, {66.2313280479,27.4339143314,45.6704667095}, {56.4895820686,23.3987510256,59.045961489}, {5.23668678881,16.1168647193,83.2935998928}, {23.4845335989,23.4845335989,78.2429125435}, {43.2352543177,43.2352543177,59.045961489}, {50.691263895,50.691263895,45.6704667095}, {64.1990066054,46.6433085545,30.4612757114}, {67.8233541008,49.2765511462,14.0290464982}, {68.7245538619,49.9313110988,-2.96645721971}, {66.8664221858,48.5812994135,-19.8428559278}, {62.323564157,45.2807198569,-35.922552248}, {48.5812994135,66.8664221858,-19.8428559278}, {49.9313110988,68.7245538619,-2.96645721971}, {49.2765511462,67.8233541008,14.0290464982}, {46.6433085545,64.1990066054,30.4612757114}, {26.2504437136,80.7905584608,-2.96645721971}, {25.9062160355,79.731134603,14.0290464982}, {24.5218384793,75.470458589,30.4612757114}, {27.4339143314,66.2313280479,45.6704667095}, {23.3987510256,56.4895820686,59.045961489}, {24.0722650693,41.6943861533,70.0507260329}, {-0.0,33.2121459216,78.2429125435}, {-0.0,48.1445301386,70.0507260329}, {-0.0,61.1438830288,59.045961489}, {-0.0,71.6882728941,45.6704667095}, {-0.0,79.3543362523,30.4612757114}, {-0.0,83.8342761307,14.0290464982}, {-0.0,84.9482202966,-2.96645721971}, {-26.2504437136,80.7905584608,-2.96645721971}, {-25.9062160355,79.731134603,14.0290464982}, {-24.5218384793,75.470458589,30.4612757114}, {-27.4339143314,66.2313280479,45.6704667095}, {-23.3987510256,56.4895820686,59.045961489}, {-24.0722650693,41.6943861533,70.0507260329}, {-13.7098240015,9.96077018859,83.2935998928}, {-23.4845335989,23.4845335989,78.2429125435}, {-43.2352543177,43.2352543177,59.045961489}, {-50.691263895,50.691263895,45.6704667095}, {-46.6433085545,64.1990066054,30.4612757114}, {-49.2765511462,67.8233541008,14.0290464982}, {-49.9313110988,68.7245538619,-2.96645721971}, {-68.7245538619,49.9313110988,-2.96645721971}, {-67.8233541008,49.2765511462,14.0290464982}, {-64.1990066054,46.6433085545,30.4612757114}, {-41.6943861533,24.0722650693,70.0507260329}, {-56.4895820686,23.3987510256,59.045961489}, {-66.2313280479,27.4339143314,45.6704667095}, {-75.470458589,24.5218384793,30.4612757114}, {-79.731134603,25.9062160355,14.0290464982}, {-80.7905584608,26.2504437136,-2.96645721971}, {-84.9482202966,5.20157830398e-15,-2.96645721971}, {-83.8342761307,5.13336889611e-15,14.0290464982}, {-79.3543362523,4.85905169449e-15,30.4612757114}, {-71.6882728941,4.38964069681e-15,45.6704667095}, {-61.1438830288,3.74398303193e-15,59.045961489}, {-48.1445301386,2.94800223654e-15,70.0507260329}, {-33.2121459216,2.03365740978e-15,78.2429125435}, {-41.6943861533,-24.0722650693,70.0507260329}, {-56.4895820686,-23.3987510256,59.045961489}, {-66.2313280479,-27.4339143314,45.6704667095}, {-75.470458589,-24.5218384793,30.4612757114}, {-79.731134603,-25.9062160355,14.0290464982}, {-80.7905584608,-26.2504437136,-2.96645721971}, {-68.7245538619,-49.9313110988,-2.96645721971}, {-67.8233541008,-49.2765511462,14.0290464982}, {-64.1990066054,-46.6433085545,30.4612757114}, {-13.7098240015,-9.96077018859,83.2935998928}, {-23.4845335989,-23.4845335989,78.2429125435}, {-43.2352543177,-43.2352543177,59.045961489}, {-50.691263895,-50.691263895,45.6704667095}, {-46.6433085545,-64.1990066054,30.4612757114}, {-49.2765511462,-67.8233541008,14.0290464982}, {-49.9313110988,-68.7245538619,-2.96645721971}, {-26.2504437136,-80.7905584608,-2.96645721971}, {-25.9062160355,-79.731134603,14.0290464982}, {-24.5218384793,-75.470458589,30.4612757114}, {-27.4339143314,-66.2313280479,45.6704667095}, {-23.3987510256,-56.4895820686,59.045961489}, {-24.0722650693,-41.6943861533,70.0507260329}, {0.0,-33.2121459216,78.2429125435}, {5.23668678881,-16.1168647193,83.2935998928}, {23.4845335989,-23.4845335989,78.2429125435}, {24.0722650693,-41.6943861533,70.0507260329}, {0.0,-48.1445301386,70.0507260329}, {0.0,-61.1438830288,59.045961489}, {0.0,-71.6882728941,45.6704667095}, {0.0,-79.3543362523,30.4612757114}, {0.0,-83.8342761307,14.0290464982}, {0.0,-84.9482202966,-2.96645721971}, {26.2504437136,-80.7905584608,-2.96645721971}, {25.9062160355,-79.731134603,14.0290464982}, {24.5218384793,-75.470458589,30.4612757114}, {27.4339143314,-66.2313280479,45.6704667095}, {23.3987510256,-56.4895820686,59.045961489}, {46.6433085545,-64.1990066054,30.4612757114}, {49.2765511462,-67.8233541008,14.0290464982}, {49.9313110988,-68.7245538619,-2.96645721971}, {48.5812994135,-66.8664221858,-19.8428559278}};

std::vector<double> biosemi_coords_to_mni(std::vector<double> v){
	arma::mat t_matrix;
	t_matrix << -4 << 0<< 0<<  84 << arma::endr
		 << 0 << 4<< 0<<  -116 << arma::endr
		 << 0 << 0<< 4<<  -56 << arma::endr
		 << 0 << 0<< 0<<  1 << arma::endr;

	std::vector<double> v1(v.begin(), v.end());
	v1.push_back(1);
	arma::rowvec o_matrix = ((t_matrix * (arma::rowvec(v1)).t()).t());
	std::vector<double> o_vec = arma::conv_to<std::vector<double>>::from(o_matrix);
	std::vector<double> o_vec1(o_vec.begin(), o_vec.end()-1);

	return o_vec1;
}

arma::mat row_norm(arma::mat& m){
	arma::mat nm(m.n_rows,1);
	//std::cout << "l2 normalization" << std::endl;
	for (int k=0; k<m.n_rows; ++k){
		float normalize_sum_squared = 0;
		for(int j=0; j<m.n_cols; ++j)
			normalize_sum_squared += std::pow(m(k,j),2);
		
		float norm = std::sqrt(normalize_sum_squared);
		nm(k,0) = norm;
	}
	return nm;
}

std::vector<arma::mat> dlegpoly(int&n, arma::mat& x){
	arma::mat nx = x.t();
	//one dimensional evaluation of the first n legendre polynomials at x

	arma::mat p(n,nx.n_cols);
	arma::mat dp1(n,nx.n_cols);
	std::vector<double> t_x;
	for(int i=0; i<nx.n_cols;++i)
		t_x.push_back((1.5 *std::pow(nx(0,i), 2) - .5));

	p.row(0) = nx;
	p.row(1) = arma::rowvec(t_x);
	dp1.row(0) = arma::ones(1,nx.n_cols);

	for(int i=2; i< n; ++i)
		p.row(i) = ((((2*i) - 1) * nx % p.row(i-1)) - ((i-1) * p.row(i-2))) / (double)i;

	for(int i=1; i< n; ++i)
		dp1.row(i) = (nx % dp1.row(i-1)) + (i*p.row(i-1));

	std::vector<arma::mat> dlegpoly{p, dp1};
	return dlegpoly;
}

arma::mat forward_model_gain_matrix_solution(std::vector<double>& rq){

	/* single dipole implementation */

	int num_dims = 3;
	int num_layers = sigmas.size();
	int num_sensors = sensor_map.size();
	arma::mat arma_re(num_sensors, num_dims);// 128 by 3
	arma::mat gain_matrix = arma::zeros(num_sensors, num_layers);
	std::vector<double> sigma_ratio(num_layers-1, 0);
	std::vector<double> sigma_ratio_sm1(num_layers-1, 0);
	double out_sphere_denom = 4.0 * PI * sigmas[num_layers-1];
	int nmax = 128; //Legendre Series Terms

	/*
		NL = length(R);                          % # of concentric sphere layers
		P = size(Rq,1);
		M = size(Re,1);
		L is 3 x nL, each column a source location
		Rq = L'; // nL x 3 ( 1 by 3)
		Re_mag = R(NL);         %(PxM)
	    Re_mag_sq = repmat(R(NL)*R(NL),P,M);  %(PxM)
	    Rq_mag = rownorm(Rq);                 %(Px1)
	    Rq_mag_sq = Rq_mag.*Rq_mag;           %(Px1)
	    Re_dot_Rq = Rq*Re';                   %(PxM)
	*/



	arma::mat arma_rq = arma::rowvec(rq);
	arma::mat arma_rq_mag = row_norm(arma_rq); //(Px1) (1 by 1)

	double re_mag = radii[num_layers-1];

	for(int i=0; i< num_sensors; ++i)
		arma_re.row(i) =  arma::rowvec(biosemi_coords_to_mni(sensor_map[i]));

	arma::mat arma_re_n = arma_re/arma::max(arma::max(arma_re));

	arma::mat arma_rq_mag_mat = arma::repmat(arma_rq_mag, 1,3);
	arma::mat arma_rq_n = arma_rq / arma_rq_mag_mat;

	for(int i=0; i<num_layers-1;++i){
		sigma_ratio[i] = sigmas[i]/sigmas[i+1]; //s = sigma_ratio
		sigma_ratio_sm1[i] = sigma_ratio[i] -1;
	}
		
	arma::mat radius_ratio = re_mag / ((arma::mat)arma::rowvec(radii)); //radius ratio
	arma::mat inv_radius_ratio = ((arma::mat)arma::rowvec(radii)) / re_mag;

	std::vector<double> twonp1(nmax,0);
	std::vector<double> n_vec(nmax,0);
	std::vector<double> f_n(nmax,0);
	for(int i=1; i<nmax+1;++i){
		twonp1[i-1] = (2*i) + 1;
		n_vec[i-1] = i;
	}

	arma::mat arma_twonp1 = arma::colvec(twonp1);	

	for(int i=0; i<nmax; ++i){
		double np1 = i + 1;
		arma::mat mc(2,2, arma::fill::eye);
		for(int j=1; j<num_layers-1; j++){
			arma::mat temp_rhs;
			temp_rhs << i+ (np1*sigma_ratio[j]) << np1 * sigma_ratio_sm1[j] * std::pow(radius_ratio[j], twonp1[i]) << arma::endr
			<< (i * sigma_ratio_sm1[j] * std::pow(inv_radius_ratio[j], twonp1[i])) << (np1 + (i*sigma_ratio[j])) << arma::endr;
			mc *= temp_rhs;
		}
		arma::mat temp_lhs;
		temp_lhs << (i * sigma_ratio_sm1[0] * std::pow(inv_radius_ratio[0], twonp1[i])) << (np1 + (i*sigma_ratio[0])) << arma::endr;
		mc.row(1) = temp_lhs * mc;
		mc /= std::pow(twonp1[i], num_layers-1);

		double denom_f_n = ((i * mc(1,1)) + (np1 * mc(1,0)));
		if(denom_f_n == 0)
			f_n[i] = 0;
		else
			f_n[i] = i/ denom_f_n;
	}
	arma::mat arma_f_n = arma::colvec(f_n);

	arma::mat onevec = arma::ones(num_sensors, 1);
	arma::mat arma_n_vec_t = (arma::rowvec(n_vec)).t();

	arma::mat arma_wtemp = ((arma_twonp1 / arma_n_vec_t) % arma_f_n) / (out_sphere_denom * std::pow(radii[num_layers-1], 2));
	
	arma::mat arma_n_vec_t_m_1 = arma_n_vec_t - 1;
	arma::mat ratio = arma_n_vec_t_m_1;

	arma_rq_n.resize(num_sensors,num_layers);
	arma_rq_n.col(0) = arma_rq_n(0,0) * onevec;
	arma_rq_n.col(1) = arma_rq_n(0,1) * onevec; 
	arma_rq_n.col(2) = arma_rq_n(0,2) * onevec;

	arma::mat cosgamma = (arma::sum(((arma_rq_n % arma_re_n).t()))).t();

	arma_rq_n = arma_rq_n.row(0);
	
	std::vector<arma::mat > pl_dp = dlegpoly(nmax, cosgamma);
	arma::mat arma_pl = pl_dp[0];
	arma::mat arma_dp = pl_dp[1];
	
	
	double ratio_base =  ((arma::mat)(arma_rq_mag / re_mag))(0,0);
	for(int i=0;i<nmax;++i)
		ratio(i,0) = std::pow(ratio_base, arma_n_vec_t_m_1(i,0));

	
	arma::mat arma_z = arma_re_n - (cosgamma * arma_rq_n);
	arma::mat arma_w = arma_wtemp % ratio;
	
	arma::mat arma_gterm1 = (arma_pl.t()) * (arma_w % arma_n_vec_t);
	
	arma::mat arma_gterm2 = (arma_dp.t()) * arma_w;
	
	arma::mat gain_matrix_right(num_sensors,num_dims);
	gain_matrix_right.col(0) = (arma_z.col(0) % arma_gterm2);
	gain_matrix_right.col(1) = (arma_z.col(1) % arma_gterm2);
	gain_matrix_right.col(2) = (arma_z.col(2) % arma_gterm2);


	gain_matrix += (arma_gterm1 * arma_rq_n) + gain_matrix_right;

    return gain_matrix;
}

arma::mat fast_forward_model_gain_matrix_solution(std::vector<double>& rq){

	int num_layers = sigmas.size();
	int num_sensors = sensor_map.size();
	arma::mat arma_re(num_sensors, num_layers);// 128 by 3
	arma::mat gain_matrix = arma::zeros(num_sensors, num_layers);

	/*
		NL = length(R);                          % # of concentric sphere layers
		P = size(Rq,1);
		M = size(Re,1);
		L is 3 x nL, each column a source location
		Rq = L'; // nL x 3 ( 1 by 3)
		Re_mag = repmat(R(NL),P,M);           %(PxM)
	    Re_mag_sq = repmat(R(NL)*R(NL),P,M);  %(PxM)
	    Rq_mag = rownorm(Rq);                 %(Px1)
	    Rq_mag_sq = Rq_mag.*Rq_mag;           %(Px1)
	    Re_dot_Rq = Rq*Re';                   %(PxM)
	*/
	
	for(int i=0; i< num_sensors; ++i)
		arma_re.row(i) =  arma::rowvec(biosemi_coords_to_mni(sensor_map[i]));
		
	std::vector<double> re_mag(num_sensors, radii[num_layers-1]);
	  
	arma::mat arma_re_mag = arma::rowvec(re_mag);


	std::vector<double> re_mag_sq(num_sensors, (radii[num_layers-1]*radii[num_layers-1])); // PxM) (3 by num_sensors)
	arma::mat arma_re_mag_sq = arma::rowvec(re_mag_sq);
	
	arma::mat arma_rq = arma::rowvec(rq);

	arma::mat arma_rq_mag = row_norm(arma_rq); //(Px1) (1 by 1)
	arma::mat arma_rq_mag_sq = arma_rq_mag % arma_rq_mag; //(Px1) (1 by 1)

	arma::mat arma_re_t = arma_re.t();
	arma::mat arma_re_dot_rq = arma_rq * arma_re_t;// (PxM) (1 by num_sensors)

    
    for(int i=0; i<num_layers;++i){
    	arma::mat arma_rq1 = mu_berg[i] * arma_rq;  //(Px3)  (1 by 3)
    	arma::mat arma_rq1_mag = mu_berg[i] * arma_rq_mag;  // (Px1) (1 by 1)
    	arma::mat arma_rq1_mag_sq = (mu_berg[i] * mu_berg[i]) * arma_rq_mag_sq; // (Px1) (1 by 1)
    	arma::mat arma_re_dot_rq1  = mu_berg[i] * arma_re_dot_rq;// (PxM) (1 by num_sensors)
    	
    	double out_sphere_denom = 4.0 * PI * sigmas[num_layers-1];
    	double inv_out_sphere_denom_lam = out_sphere_denom/lam_berg[i];
    	arma::mat term = 1/(inv_out_sphere_denom_lam * arma_rq1_mag_sq); // (Px1) (1 by 1)


    	std::vector<int> nx;
    	for(int i=0;i<arma_rq1_mag.n_rows;++i){
    		if(arma_rq1_mag(i,0) > radii[num_layers-1])
    			nx.push_back(i);
    	}
    	if(nx.size()> 0){
    		std::cout << "berg dipole outside of sphere" << std::endl;
    		for(int i =0; i< nx.size(); ++i){
    			arma::mat rq1_temp = arma_rq1.row(nx[i]);
    			arma::mat rq1_temp_mag_sq = row_norm(rq1_temp) % row_norm(rq1_temp);
    			arma_rq1.row(nx[i]) = (radii[num_layers-1] * radii[num_layers-1])/arma::repmat(rq1_temp_mag_sq, 1, 3);
    			arma::mat sub_rq1 = arma_rq1.row(nx[i]);
    			arma_rq1_mag.row(nx[i]) = row_norm(sub_rq1);
    			arma_rq1_mag_sq(nx[i],0) = arma_rq1_mag(nx[i],0) * arma_rq1_mag(nx[i],0);
    			arma_re_dot_rq1.row(nx[i]) = ((radii[num_layers-1] * radii[num_layers-1]) * arma_re_dot_rq1.row((nx[i]))) / arma::repmat(rq1_temp_mag_sq, 1, num_sensors);
    			
    			arma::mat normed_rq1_temp = (radii[num_layers-1]/row_norm(rq1_temp));
    			term.row(nx[i]) = normed_rq1_temp.t() % term.row(nx[i]);
    		}
    	}

    	//Calculation of Forward Gain Matrix Contribution due to K-th Berg Dipole

    	arma_rq1_mag = arma::repmat(arma_rq1_mag, 1, num_sensors); // (PxM) (1 by num_sensors)
    	arma_rq1_mag_sq = arma::repmat(arma_rq1_mag_sq, 1, num_sensors); // (PxM) (1 by num_sensors)
    	term = arma::repmat(term, 1, num_sensors); // (PxM) (1 by num_sensors)

    	/* D_MAG AND D_MAG^3 CALC BELOW */
    	arma::mat re_temp_kern_trans = arma::repmat(arma_re, 1, 1).t();
    	arma::mat re_temp_kern_trans_rn = arma::reshape(re_temp_kern_trans, 3, num_sensors).t();
    	arma::mat temp_norm_re_temp_rq1  = re_temp_kern_trans_rn - arma::repmat(arma_rq1, num_sensors, 1);
    	arma::mat d_mag = arma::reshape(row_norm(temp_norm_re_temp_rq1), 1, num_sensors); // (PxM) (1 by num_sensors)
    	arma::mat d_mag_cub = d_mag % d_mag % d_mag; // (PxM) (1 by num_sensors)
    	/* D_MAG AND D_MAG^3 CALC ABOVE */

    	arma::mat f_scalar = d_mag % ( ((arma_re_mag % d_mag) + arma_re_mag_sq) - arma_re_dot_rq1 ); // (PxM) (1 by num_sensors)

    	arma::mat c1 = term % (((2 *(arma_re_dot_rq1 - arma_rq1_mag_sq) / d_mag_cub )) + ((1/d_mag) - (1/arma_re_mag))); // (PxM) (1 by num_sensors)
    	
    	arma::mat c2 = term % ((2 / d_mag_cub) + ((d_mag + arma_re_mag) / (arma_re_mag % f_scalar))); // (PxM) (1 by num_sensors)

    	gain_matrix += (arma::reshape(arma::repmat((c1 - (c2 % arma_re_dot_rq1)).t(), 3, 1), num_sensors, 3) % arma::repmat(arma::reshape(arma_rq1.t(), 1, 3), num_sensors, 1)) + (arma::reshape(arma::repmat((c2 % arma_rq1_mag_sq).t(), 3, 1), num_sensors, 3) % arma::repmat(arma_re, 1, 1));
    }

    return gain_matrix;
}
int main(int argc, char *argv[]){

	std::vector<double> rq {-6,-60,18};
	arma::mat lead_field_potential = forward_model_gain_matrix_solution(rq);
	arma::mat fs = lead_field_potential.t();
	for(int i=0; i< fs.n_cols;++i){
		double mag_l = std::sqrt(std::pow(fs(0,i),2) + std::pow(fs(1,i),2) + std::pow(fs(2,i),2));
		fs.col(i) = fs.col(i)/mag_l;
	}
	std::cout << fs << std::endl;
	return 0;
}