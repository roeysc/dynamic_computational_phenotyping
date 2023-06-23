

functions {

    real partial_sum(int[,,] choice_itc, int start, int end, int[,,] choice_nc, real[,] weber_nc, int[,] idx_nc_obs, 
                    int[,,] deltaM, real[,,] TotalS,  int[,] Tr_nc, int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {

 	        if (idx_nc_obs[s,w] != 0) { 	 
 	            vector[Tr_nc[s,w]] z = - to_vector(deltaM[s,w,:Tr_nc[s,w]]) ./ (weber_nc[s,w] * to_vector(TotalS[s,w,:Tr_nc[s,w]]));
	            lt += bernoulli_lpmf(choice_nc[s,w,:Tr_nc[s,w]] | Phi(z));	
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
    int<upper=1> choice_itc[N,W,Tr_max_itc];     // choice itc - 0 for instant reward, 1 for delayed reward - subj x weeks x trials

    // Numerocity comparison 
    int P_nc;
    int<lower=1> W_nc_obs[N];                                          // Number of weeks WITH data
    int<lower=0> idx_nc_obs[N,W];                              // Indices of weeks WITH data
    int<lower=0> Tr_max_nc;                                          // Max number of trials across subjects across weeks
    int<lower=0> Tr_nc[N, W];                                // Number of trials for each subj for each week
    int deltaM[N, W, Tr_max_nc];                              // Mu of dots - subj x weeks x trials
    real TotalS[N, W, Tr_max_nc];                             // Sd of dots - subj x weeks x trials
    int choice_nc[N, W, Tr_max_nc];                           // choice nc - subj x weeks x trials

}

transformed data {
    int num_par =P_nc;
}

parameters {
    real mu_pr;
    real<lower=0> sigma_pr;
    real mu_pr_sub[N];
    real<lower=0> sigma_pr_r;

    real weber_nc_pr[N,W];                  // nc

}

transformed parameters {
    real<lower=0> weber_nc[N,W];

    for (n in 1:N) {
        for (w in 1:W) {
            weber_nc[n,w] = exp(mu_pr_sub[n] + sigma_pr_r*weber_nc_pr[n,w]);          // nc
        }
    }

}

model {
    // hyper parameters
    mu_pr  ~ normal(0, 1);
    sigma_pr ~ normal(0, 1);

    // subject level parameters
    to_vector(mu_pr_sub) ~ normal(mu_pr, sigma_pr);

    sigma_pr_r ~ normal(0, 1);

    // individual parameters w/ Matt trick
    to_vector(to_matrix(weber_nc_pr)) ~ normal(0, 1);

    	target += reduce_sum(partial_sum, choice_itc, 1, choice_nc, weber_nc, idx_nc_obs, deltaM, TotalS, Tr_nc, W);
}    

 generated quantities {
    // For posterior predictive check
    real y_pred_all_weeks_nc[Tr_max_nc, W,N];

    real log_lik_all_subs[N, W];

    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {
            for (t in 1:Tr_max_nc) {
                y_pred_all_weeks_nc[t,w,n] = -2;
            }
        }
    }

    // here start predictive checks and logliklhd calculations
    for (n in 1:N) {
        for (w in 1:W) {

            log_lik_all_subs[n,w] = 0;   

            if (idx_nc_obs[n,w] != 0) {
                for (t in 1:Tr_nc[n,w]) {
                 	real z = - (deltaM[n,w,t]) ./ (weber_nc[n,w] * TotalS[n,w,t]);
                    // generate posterior prediction for current trial
                    y_pred_all_weeks_nc[t,w,n] = bernoulli_rng(normal_cdf(0, deltaM[n,w, t], weber_nc[n,w] * TotalS[n,w, t]));   
                    log_lik_all_subs[n,w] +=  bernoulli_lpmf(choice_nc[n,w,t] | Phi(z));    
                }
            }                      
        }
    }
 }


