

functions {

    real partial_sum(int[,,] choice_itc, int start, int end,

			        int[,,] choice_cd, real[,] criterion_cd, real[,] sigma_cd, real[,] criterion_slope_cd, real[,] sigma_slope_cd, 
			        int[,] idx_cd_obs, int[,,] Nb_cd, int[,,] Tar_cd,  real[,,] delta_cd, int[,]  Tr_cd, int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {

 	        if (idx_cd_obs[s,w] != 0) { 
                for (t in 1:Tr_cd[s,w]) {
                    // Figure out what block this is based on the number of targets
        	        real factor = Nb_cd[s,w,t]-3; // // Factor is the factor by which we need to multiply sigma and critrion, which are assumed to increase linearly with the number of items. The baseline is 3, which is the minimal number of items in the task.

                    real f = 0.5 * (erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor))/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))) + erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor))/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))));       

                    // generate posterior prediction for current trial
                    if (delta_cd[s,w,t]==0){
                	    real h = 0.5 * (erfc((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))) + erfc((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))));
                        real p = pow((1-h), Tar_cd[s,w,t]) * pow((1-f),(Nb_cd[s,w,t]-Tar_cd[s,w,t]));   
                        lt += bernoulli_lpmf(choice_cd[s,w,t] | 1.0-p); // Becuase choice_cd is what's being sliced over, it's indexed with n and not s
                	} else {
                	    real h = 0.5 * (erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)-1)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))) + erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)+1)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))));
                        real p = pow((1-h), Tar_cd[s,w,t]) * pow((1-f),(Nb_cd[s,w,t]-Tar_cd[s,w,t]));   
                        lt += bernoulli_lpmf(choice_cd[s,w,t] | 1.0-p);
                    }
                }
	        }	         	    	
		}
	}
	return lt;
	}
}


data {

    // Intertemporal choice
    int N;                                                          // Number of subjects
    int W;                                                          // Number of weeks (typically 12)
    
    int P_itc;
    int<lower=0> idx_itc_obs[N,W];                            // Indices of weeks WITH data
    int<lower=0> Tr_max_itc;                           // Max number of trials across subjects across weeks
    int<lower=0> Tr_itc[N, W];                              // Number of trials for each subj for each week
    real<lower=0> amount_later[N,W,Tr_max_itc];            // Amount later - subj x weeks x trials      
    real<lower=0> amount_sooner[N,W,Tr_max_itc];           // Amount sooner - subj x weeks x trials
    real<lower=0> delay_later[N,W,Tr_max_itc];             // Delay later - subj x weeks x trials
    int<lower=-1, upper=1> choice_itc[N,W,Tr_max_itc];     // choice itc - 0 for instant reward, 1 for delayed reward - subj x weeks x trials

    // Change detection 
    int P_cd;
    int<lower=0> idx_cd_obs[N,W];                              // Indices of weeks WITH data
    int<lower=0> Tr_max_cd;                                          // Max number of trials across subjects across weeks
    int<lower=0> Tr_cd[N, W];                                // Number of trials for each subj for each week
    int<lower=0> Nb_cd[N, W, Tr_max_cd];                              // Number of objects inper trial
    int<lower=0> Tar_cd[N, W, Tr_max_cd];                              // Number of targets per trial
    real<lower=0> delta_cd[N, W, Tr_max_cd];                             // Whether or not a change occured per trial
    int choice_cd[N, W, Tr_max_cd];                           // choice cd - subj x weeks x trials 
}

transformed data {
    int num_par = P_cd;
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];

    real criterion_cd_pr[N,W]; //    4 criteria, one for each block
    real sigma_cd_pr[N,W];     // 4 sigmas, one for each block
    real criterion_slope_cd_pr[N,W]; //    4 criteria, one for each block
    real sigma_slope_cd_pr[N,W];     // 4 sigmas, one for each block

}

transformed parameters {
    real criterion_cd[N,W]; //    4 criteria, one for each block
    real sigma_cd[N,W];     // 4 sigmas, one for each block
    real criterion_slope_cd[N,W]; //    4 criteria, one for each block
    real sigma_slope_cd[N,W];     // 4 sigmas, one for each block

    for (n in 1:N) {
        for (w in 1:W) {
            criterion_cd[n,w] = exp(mu_pr_sub[1,n] + sigma_pr_r[1]*criterion_cd_pr[n,w]);
            criterion_slope_cd[n,w] = exp(mu_pr_sub[2,n] + sigma_pr_r[2]*criterion_slope_cd_pr[n,w]);
            sigma_cd[n,w] = exp(mu_pr_sub[3,n] + sigma_pr_r[3]*sigma_cd_pr[n,w]);
            sigma_slope_cd[n,w] = exp(mu_pr_sub[4,n] + sigma_pr_r[4]*sigma_slope_cd_pr[n,w]);
        }
    }
}

model {
    // hyper parameters
    mu_pr  ~ normal(0, 1);
    sigma_pr ~ normal(0, 1);

    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr[p], sigma_pr[p]);
    }

    to_vector(sigma_pr_r) ~ normal(0, 1);

    // individual parameters w/ Matt trick
    to_vector(to_matrix(criterion_cd_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(criterion_slope_cd_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(sigma_cd_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(sigma_slope_cd_pr)) ~ normal(0, 1.0);

    	target += reduce_sum(partial_sum, choice_itc, 1, choice_cd, criterion_cd, sigma_cd, criterion_slope_cd, 
    	                    sigma_slope_cd, idx_cd_obs, Nb_cd, Tar_cd, delta_cd, Tr_cd, W);
}    

 generated quantities {
    // For posterior predictive check
    real y_pred_all_weeks_cd[Tr_max_cd,W,N];

    real log_lik_all_subs[N, W];

    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {         
            for (t in 1:Tr_max_cd) {
                y_pred_all_weeks_cd[t,w,n] = -2;
            }

        }
    }

    // here start predictive checks and logliklhd calculations
    for (n in 1:N) {
        for (w in 1:W) {

            log_lik_all_subs[n,w] = 0;

            if (idx_cd_obs[n,w] != 0) {
                for (t in 1:Tr_cd[n,w]) {
                    // Skip trials with more than 1 target
                    //if (Tar_cd[n,w,t]>1) continue;  // I decided against skipping the trials with more than one targte, because the model can take care of that

    	            real factor = Nb_cd[n,w,t]-3; // // Factor is the factor by which we need to multiply sigma and critrion, which are assumed to increase linearly with the number of items. The baseline is 3, which is the minimal number of items in the task.

                    real f = 0.5 * (erfc(((criterion_cd[n,w]+criterion_slope_cd[n,w]*factor))/(sqrt(2)*(sigma_cd[n,w]+sigma_slope_cd[n,w]*factor))) + erfc(((criterion_cd[n,w]+criterion_slope_cd[n,w]*factor))/(sqrt(2)*(sigma_cd[n,w]+sigma_slope_cd[n,w]*factor))));       

                    // generate posterior prediction for current trial
                    if (delta_cd[n,w,t]==0){
            	        real h = 0.5 * (erfc((criterion_cd[n,w]+criterion_slope_cd[n,w]*factor)/(sqrt(2)*(sigma_cd[n,w]+sigma_slope_cd[n,w]*factor))) + erfc((criterion_cd[n,w]+criterion_slope_cd[n,w]*factor)/(sqrt(2)*(sigma_cd[n,w]+sigma_slope_cd[n,w]*factor))));
                        real p = pow((1-h), Tar_cd[n,w,t]) * pow((1-f),(Nb_cd[n,w,t]-Tar_cd[n,w,t]));   
                        y_pred_all_weeks_cd[t,w,n] = bernoulli_rng(1-p);  
                        log_lik_all_subs[n,w] += bernoulli_lpmf(choice_cd[n,w,t] | 1-p);
            	    } else {
            	        real h = 0.5 * (erfc(((criterion_cd[n,w]+criterion_slope_cd[n,w]*factor)-1)/(sqrt(2)*(sigma_cd[n,w]+sigma_slope_cd[n,w]*factor))) + erfc(((criterion_cd[n,w]+criterion_slope_cd[n,w]*factor)+1)/(sqrt(2)*(sigma_cd[n,w]+sigma_slope_cd[n,w]*factor))));
                        real p = pow((1-h), Tar_cd[n,w,t]) * pow((1-f),(Nb_cd[n,w,t]-Tar_cd[n,w,t]));   
                        y_pred_all_weeks_cd[t,w,n] = bernoulli_rng(1-p);  
                        log_lik_all_subs[n,w] += bernoulli_lpmf(choice_cd[n,w,t] | 1-p);
                   }
                }
            }                     
        }
    }
 }


