

functions {

    real partial_sum(int[,,] choice_itc, int start, int end, real[,,] delay_later, real[,,] amount_later, real[,,] amount_sooner, 
                    real[,] itc_k, real[,] itc_beta, int[,] idx_itc_obs, int[,] Tr_itc, int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {

		    if (idx_itc_obs[s,w] != 0) {

                vector[Tr_itc[s,w]] ev_later   = to_vector(amount_later[s,w,:Tr_itc[s,w]])  ./ (1 + itc_k[s,w] * to_vector(delay_later[s,w,:Tr_itc[s,w]]));
 	            vector[Tr_itc[s,w]] ev_sooner  = to_vector(amount_sooner[s,w,:Tr_itc[s,w]]);
      			lt += bernoulli_logit_lpmf(choice_itc[n,w,:Tr_itc[s,w]] | itc_beta[s,w] * (ev_later - ev_sooner));
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
}

transformed data {
    int num_par = P_itc;
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];

    real itc_k_pr[N,W];                     // itc
    real itc_beta_pr[N,W];                  // itc

}

transformed parameters {
    real<lower=0> itc_k[N,W];           // itc
    real<lower=0> itc_beta[N,W];        // itc
    for (n in 1:N) {
        for (w in 1:W) {
            itc_k[n,w] = exp(mu_pr_sub[1,n] + sigma_pr_r[1]*itc_k_pr[n,w]);                // itc
            itc_beta[n,w] = exp(mu_pr_sub[2,n] + sigma_pr_r[2]*itc_beta_pr[n,w]);          // itc
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
    to_vector(to_matrix(itc_k_pr)) ~ normal(0, 1);
    to_vector(to_matrix(itc_beta_pr)) ~ normal(0, 1);

    	target += reduce_sum(partial_sum, choice_itc, 1, delay_later, amount_later, amount_sooner, itc_k, itc_beta, 
    	                    idx_itc_obs, Tr_itc, W);
}    


 generated quantities {
    // For posterior predictive check
    real y_pred_all_weeks_itc[Tr_max_itc, W,N];

    real log_lik_all_subs[N, W];

    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {
            for (t in 1:Tr_max_itc) {
                y_pred_all_weeks_itc[t, w,n] = -2;
            }
        }
    }

    // here start predictive checks and logliklhd calculations
    for (n in 1:N) {
        for (w in 1:W) {

            log_lik_all_subs[n,w] = 0;

            if (idx_itc_obs[n,w] != 0) {
                for (t in 1:(Tr_itc[n,w])) {
                    real  ev_later  = amount_later[n,w,t]  / (1 + itc_k[n,w] * delay_later[n,w,t]); # Once upon a time we also divided k_itc by 7 here, but we no longer do this, to avoid confusion.
                    real ev_sooner  = amount_sooner[n,w,t];
                    // generate posterior prediction for current trial (and repeat it for 12 "weeks" to have 12 simulated datasets for each session)
                    y_pred_all_weeks_itc[t,w,n] = bernoulli_rng(inv_logit(itc_beta[n,w] * (ev_later - ev_sooner)));
                    log_lik_all_subs[n,w] += bernoulli_logit_lpmf(choice_itc[n,w,t] | itc_beta[n,w] * (ev_later - ev_sooner));
                }
            }                   
        }
    }
 }


