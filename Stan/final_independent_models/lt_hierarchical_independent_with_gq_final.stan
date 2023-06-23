

functions {

    real partial_sum(int[,,] choice_itc, int start, int end,  int[,,] choice_lt, real[,,] hi_p_lt, real[,,] hi_narr_lt,
                    real[,,] lo_narr_lt, real[,,] hi_wide_lt, real[,,] lo_wide_lt, int[,] Tr_lt, real[,] risk_lt,
                    real[,] beta_lt, int[,] idx_lt_obs, int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {

           if (idx_lt_obs[s,w] != 0) { 	 

        		for (t in 1:Tr_lt[s,w]) {
          			real evNarr;    // evNarr, evWide, pSafe can be a scalar to save memory and increase speed.
          			real evWide;  // they are left as arrays as an example for RL models.
        			real evNarr_var;    // varaiance of evNarr
          			real evWide_var;  // varaiance of evWode
          			real pSafe;

          			evNarr  = (1-hi_p_lt[s,w,t])*pow(lo_narr_lt[s,w,t], risk_lt[s,w]) + hi_p_lt[s,w,t]*pow(hi_narr_lt[s,w,t], risk_lt[s,w]);
          			evWide  = (1-hi_p_lt[s,w,t])*pow(lo_wide_lt[s,w,t], risk_lt[s,w]) + hi_p_lt[s,w,t]*pow(hi_wide_lt[s,w,t], risk_lt[s,w]);
                    evNarr_var = (1-hi_p_lt[s,w,t])*pow(pow(lo_narr_lt[s,w,t], risk_lt[s,w])-evNarr, 2) + hi_p_lt[s,w,t]*pow(pow(hi_narr_lt[s,w,t], risk_lt[s,w])-evNarr, 2);
                    evWide_var = (1-hi_p_lt[s,w,t])*pow(pow(lo_wide_lt[s,w,t], risk_lt[s,w])-evWide, 2) + hi_p_lt[s,w,t]*pow(pow(hi_wide_lt[s,w,t], risk_lt[s,w])-evWide, 2);

          			pSafe  = inv_logit(beta_lt[s,w] * (evNarr - evWide)/sqrt(evNarr_var+evWide_var));

          			lt += bernoulli_lpmf(choice_lt[s,w,t] | pSafe); // Here we haveindex with n and not with s because choice_lt is the input argument that's being sliced over
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
    int<upper=1> choice_itc[N,W,Tr_max_itc];     // choice itc - 0 for instant reward, 1 for delayed reward - subj x weeks x trials

    int P_lt;
    int<lower=0> idx_lt_obs[N,W];                              // Indices of weeks WITH data
    int<lower=0> Tr_max_lt;                                          // Max number of trials across subjects across weeks
    int<lower=0> Tr_lt[N, W];                                // Number of trials for each subj for each week
    real<lower=0> hi_p_lt[N, W,Tr_max_lt];
    real<lower=0> hi_narr_lt[N, W, Tr_max_lt];
    real<lower=0> lo_narr_lt[N, W, Tr_max_lt];
    real<lower=0> hi_wide_lt[N, W, Tr_max_lt];
    real<lower=0> lo_wide_lt[N, W, Tr_max_lt];
    int<upper=1> choice_lt[N, W, Tr_max_lt]; // Risky=0; Safe=1. Roey: Why is -1 an option here?   
}

transformed data {
    int num_par = P_lt;
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];

    real risk_lt_pr[N,W];               // lt
    real beta_lt_pr[N,W];               // lt
}

transformed parameters {
    real<lower=0> risk_lt[N,W];        // lt
    real<lower=0> beta_lt[N,W];       // lt

    for (n in 1:N) {
        for (w in 1:W) {
            risk_lt[n,w] = exp(mu_pr_sub[1,n] + sigma_pr_r[1]*risk_lt_pr[n,w]);
            beta_lt[n,w] = exp(mu_pr_sub[2,n] + sigma_pr_r[2]*beta_lt_pr[n,w]);
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
    to_vector(to_matrix(risk_lt_pr)) ~ normal(0, 1);
    to_vector(to_matrix(beta_lt_pr)) ~ normal(0, 1);


    	target += reduce_sum(partial_sum, choice_itc, 1, choice_lt, hi_p_lt, hi_narr_lt, lo_narr_lt, hi_wide_lt, 
    	                    lo_wide_lt, Tr_lt, risk_lt, beta_lt, idx_lt_obs, W);
}    

 generated quantities {
    // For posterior predictive check
    real y_pred_all_weeks_lt[Tr_max_lt,W,N];

    real log_lik_all_subs[N, W];

    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {
            for (t in 1:Tr_max_lt) {
                y_pred_all_weeks_lt[t, w, n] = -2;
            }            
        }
    }

    // here start predictive checks and logliklhd calculations
    for (n in 1:N) {
        for (w in 1:W) {

            log_lik_all_subs[n,w] = 0;

            if (idx_lt_obs[n,w] != 0) {
                for (t in 1:Tr_lt[n,w]) {
                    real evNarr;    // evNarr, evWide, pSafe can be a scalar to save memory and increase speed.
                    real evWide;  // they are left as arrays as an example for RL models.
       			    real evNarr_var;    // varaiance of evNarr
      			    real evWide_var;  // varaiance of evWode
                    real pSafe;

                    evNarr  = (1-hi_p_lt[n,w,t])*pow(lo_narr_lt[n,w,t], risk_lt[n,w]) + hi_p_lt[n,w,t]*pow(hi_narr_lt[n,w,t], risk_lt[n,w]);
                    evWide  = (1-hi_p_lt[n,w,t])*pow(lo_wide_lt[n,w,t], risk_lt[n,w]) + hi_p_lt[n,w,t]*pow(hi_wide_lt[n,w,t], risk_lt[n,w]);
                    evNarr_var = (1-hi_p_lt[n,w,t])*pow(pow(lo_narr_lt[n,w,t], risk_lt[n,w])-evNarr, 2) + hi_p_lt[n,w,t]*pow(pow(hi_narr_lt[n,w,t], risk_lt[n,w])-evNarr, 2);
                    evWide_var = (1-hi_p_lt[n,w,t])*pow(pow(lo_wide_lt[n,w,t], risk_lt[n,w])-evWide, 2) + hi_p_lt[n,w,t]*pow(pow(hi_wide_lt[n,w,t], risk_lt[n,w])-evWide, 2);

        		    pSafe  = inv_logit(beta_lt[n,w] * (evNarr - evWide)/sqrt(evNarr_var+evWide_var));

                    log_lik_all_subs[n,w] += bernoulli_lpmf(choice_lt[n,w,t] | pSafe);

                    // generate posterior prediction for current trial
                    y_pred_all_weeks_lt[t,w,n] = bernoulli_rng(pSafe);
                }
            }                   
        }
    }
 }


