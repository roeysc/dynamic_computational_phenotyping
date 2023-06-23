

functions {

    real partial_sum(int[,,] choice_itc, int start, int end, int[,,] y_smb, real[,,] params_smb, int[,] idx_smb_obs, 
                        real[,,,] x_smb, int[,] Tr_smb,  int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {

	        if (idx_smb_obs[s,w] != 0) {                                    // if week exists

	            matrix[Tr_smb[s,w], 3] X_smb = to_matrix(x_smb[s,w,:Tr_smb[s,w],]);
	            lt += bernoulli_lpmf(y_smb[s,w,:Tr_smb[s,w]] | Phi(X_smb * to_vector(params_smb[s,w,])));
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

    // Two-armed bandit
    int<lower=1> P_smb; // number of predictors
    int<lower=0> idx_smb_obs[N,W];                            // Indices of weeks WITH data
    int<lower=1> Tr_max_smb; // max number of trials
    int<lower=0, upper=Tr_max_smb> Tr_smb[N, W]; 
    real x_smb[N, W, Tr_max_smb, P_smb];  // predictor matrix
    int y_smb[N, W, Tr_max_smb];  // Outcome - sub X trials

}

transformed data {
    int num_par = P_smb;
}

parameters {

    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];

    real params_smb_pr[N, W, 3];            // smb 3 params

}

transformed parameters {
    real params_smb[N, W, 3];           // smb

    for (n in 1:N) {
        for (w in 1:W) {
            params_smb[n,w,1] = (mu_pr_sub[1,n] + sigma_pr_r[1]*params_smb_pr[n,w,1]);     // smb                                  
            params_smb[n,w,2]  = (mu_pr_sub[2,n] + sigma_pr_r[2]*params_smb_pr[n,w,2]);    // smb
            params_smb[n,w,3] = (mu_pr_sub[3,n] + sigma_pr_r[3]*params_smb_pr[n,w,3]);     // smb
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
    to_vector(to_matrix(params_smb_pr[,,1])) ~ normal(0, 1);
    to_vector(to_matrix(params_smb_pr[,,2])) ~ normal(0, 1);
    to_vector(to_matrix(params_smb_pr[,,3])) ~ normal(0, 1);

    	target += reduce_sum(partial_sum, choice_itc, 1, y_smb, params_smb, idx_smb_obs, x_smb, Tr_smb, W);
}    

 generated quantities {
    // For posterior predictive check
    real y_pred_all_weeks_smb[Tr_max_smb, W, N];

    real log_lik_all_subs[N, W];

    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {
            for (t in 1:Tr_max_smb) {
                y_pred_all_weeks_smb[t, w, n] = -2;
            }
        }
    }

    // here start predictive checks and logliklhd calculations
    for (n in 1:N) {
        for (w in 1:W) {

            log_lik_all_subs[n,w] = 0;

            if (idx_smb_obs[n,w] != 0) {                                // if week exists
                    matrix[Tr_smb[n,w], P_smb] X_smb;
                    X_smb = to_matrix(x_smb[n,w,:Tr_smb[n,w],]);
                    y_pred_all_weeks_smb[:Tr_smb[n,w], w, n] = bernoulli_rng(Phi(X_smb * to_vector(params_smb[n,w,])));   
                    log_lik_all_subs[n,w] += bernoulli_lpmf(y_smb[n,w,:Tr_smb[n,w]] | Phi(X_smb * to_vector(params_smb[n,w,])));               
            }                         
        }
    }
 }


