

functions {

    real partial_sum(int[,,] choice_itc, int start, int end,  real[,,,] outcome_gng, int[,] idx_gng_obs,
                    int[,,,] pressed_gng, int[,,,] cue_gng,  int[,,] Tr_gng, int Bl, real[,] b_gng, real[,] pi_gng, 
                    real[,] xi_gng, real[,] ep_gng, real[,] rho_gng,  int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {

	        if (idx_gng_obs[s,w] != 0) {                                // if week exists

	                for (b in 1:Bl) {

	                    vector[4] wv_g  = rep_vector(0.0, 4);  // action weight for go
	                    vector[4] wv_ng = rep_vector(0.0, 4); // action weight for nogo
	                    vector[4] qv_g = rep_vector(0.0, 4);  // Q value for go
	                    vector[4] qv_ng = rep_vector(0.0, 4); // Q value for nogo
	                    vector[4] sv = rep_vector(0.0, 4);    // stimulus value
	                    vector[4] pGo = rep_vector(0.0, 4);   // prob of go (press)

	                    for (t in 1:Tr_gng[s,w,b]) {
	                        wv_g[cue_gng[s,w,b,t]]  = qv_g[cue_gng[s,w,b,t]] + b_gng[s,w] + pi_gng[s,w] * sv[cue_gng[s,w,b,t]];
	                        wv_ng[cue_gng[s,w,b,t]] = qv_ng[cue_gng[s,w,b,t]];  // qv_ng is always equal to wv_ng (regardless of action)
	                        pGo[cue_gng[s,w, b,t]]   = inv_logit(wv_g[cue_gng[s,w, b,t]] - wv_ng[cue_gng[s,w,b,t]]);
	                        pGo[cue_gng[s,w,b,t]]   *= (1 - xi_gng[s,w]);
	                        pGo[cue_gng[s,w,b,t]]   += xi_gng[s,w]/2;

	                        lt += bernoulli_lpmf(pressed_gng[s,w,b,t] | pGo[cue_gng[s,w,b,t]]);

	                        // after receiving feedback, update sv[t + 1]
	                        sv[cue_gng[s,w,b,t]] += ep_gng[s,w] * (rho_gng[s,w] * outcome_gng[s,w,b,t] - sv[cue_gng[s,w,b,t]]);
	                        // update action values
	                        if (pressed_gng[s,w,b,t] == 1) { // update go value
	                            qv_g[cue_gng[s,w,b,t]] += ep_gng[s,w] * (rho_gng[s,w] * outcome_gng[s,w,b,t] - qv_g[cue_gng[s,w,b,t]]);
                            }
	                        else if  (pressed_gng[s,w,b,t] == 0) { // update no-go value
	                            qv_ng[cue_gng[s,w,b,t]] += ep_gng[s,w] * (rho_gng[s,w] * outcome_gng[s,w,b,t] - qv_ng[cue_gng[s,w,b,t]]);
                            }
	                    } // end of t loop
	                } // end of b loop
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

    int P_gng;
    int<lower=0> idx_gng_obs[N,W];                            // Indices of weeks WITH data
    int<lower=0> Bl;
    int<lower=1> Tr_max_gng;
    int<lower=0, upper=Tr_max_gng> Tr_gng[N,W,Bl];
    int<lower=0, upper=4> cue_gng[N,W,Bl,Tr_max_gng];
    int<upper=1> pressed_gng[N,W,Bl, Tr_max_gng];
    real outcome_gng[N, W, Bl, Tr_max_gng];  

}

transformed data {
    int num_par = P_gng;
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];

    real xi_gng_pr[N,W];                 // gng, noise
    real ep_gng_pr[N,W];                 // gng, learning rate
    real rho_gng_pr[N,W];                // gng, rho
    real b_gng_pr[N,W];
    real pi_gng_pr[N,W];

}

transformed parameters {
    real b_gng[N,W];
    real<lower=0> pi_gng[N,W];
    real<lower=0, upper=1> xi_gng[N,W];                 // gng, noise
    real<lower=0, upper=1> ep_gng[N,W];                 // gng, learning rate
    real<lower=0> rho_gng[N,W];                // gng, rho

    for (n in 1:N) {
        for (w in 1:W) {
            xi_gng[n,w] = Phi(mu_pr_sub[1,n] + sigma_pr_r[1]*xi_gng_pr[n,w]);
            ep_gng[n,w] = Phi(mu_pr_sub[2,n] + sigma_pr_r[2]*ep_gng_pr[n,w]);
            b_gng[n,w] = (mu_pr_sub[3,n] + sigma_pr_r[3]*b_gng_pr[n,w]);
            pi_gng[n,w] = exp(mu_pr_sub[4,n] + sigma_pr_r[4]*pi_gng_pr[n,w]);
            rho_gng[n,w] = exp(mu_pr_sub[5,n] + sigma_pr_r[5]*rho_gng_pr[n,w]);
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
    to_vector(to_matrix(xi_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(ep_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(rho_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(b_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(pi_gng_pr)) ~ normal(0, 1);

    	target += reduce_sum(partial_sum, choice_itc, 1, outcome_gng, idx_gng_obs, pressed_gng, cue_gng, Tr_gng, Bl,
    	                    b_gng, pi_gng, xi_gng, ep_gng, rho_gng, W);
}    

 generated quantities {
    // For posterior predictive check
    real choice_ppc_gng[Tr_max_gng, Bl, W,N];
    real pGo_ppc_gng[Tr_max_gng, Bl, W,N];

    real log_lik_all_subs[N, W];

    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {
            for (t in 1:Tr_max_gng) {
                for (bl in 1:Bl) {
                    choice_ppc_gng[t,bl,w,n] = -2;
                    pGo_ppc_gng[t,bl,w,n] = -2;
                }
            }
        }
    }

    // here start predictive checks and logliklhd calculations
    for (n in 1:N) {
        for (w in 1:W) {

            log_lik_all_subs[n,w] = 0;

            if (idx_gng_obs[n,w] != 0) {                                // if week exists

                for (bl in 1:Bl) {
                // Treat each simulation (w 1-12) and block (bl 1-3) as an independent process, where we update the wv's and qv's independently from other simulations 
                    vector[4] wv_g  = rep_vector(0.0, 4);
                    vector[4] wv_ng = rep_vector(0.0, 4);
                    vector[4] qv_g  = rep_vector(0.0, 4);
                    vector[4] qv_ng = rep_vector(0.0, 4);
                    vector[4] sv    = rep_vector(0.0, 4);
                    vector[4] pGo   = rep_vector(0.0, 4);   // prob of go (press)

                    for (t in 1:Tr_gng[n,w,bl]) {
                        wv_g[cue_gng[n,w,bl,t]]  = qv_g[cue_gng[n,w,bl,t]] + b_gng[n,w] + pi_gng[n,w] * sv[cue_gng[n,w,bl, t]];
                        wv_ng[cue_gng[n,w,bl,t]] = qv_ng[cue_gng[n,w,bl, t]];  // qv_ng is always equal to wv_ng (regardless of action)
                        pGo[cue_gng[n,w,bl, t]]   = inv_logit(wv_g[cue_gng[n,w,bl, t]] - wv_ng[cue_gng[n,w,bl,t]]);
                        {  // noise
                        pGo[cue_gng[n,w,bl,t]]   *= (1 - xi_gng[n,w]);
                        pGo[cue_gng[n,w,bl,t]]   += xi_gng[n,w]/2;
                        }

                        // Keep track of the pGo value in this trial
                        pGo_ppc_gng[t,bl,w,n] = pGo[cue_gng[n,w,bl,t]];

                        // And keep track of the predicted choice made in this trial
                        choice_ppc_gng[t,bl,w,n] = bernoulli_rng(pGo[cue_gng[n,w,bl, t]]);

                        log_lik_all_subs[n,w] += bernoulli_lpmf(pressed_gng[n,w, bl,t] | pGo[cue_gng[n,w,bl, t]]);

                        // Update values based on  the action selected by the real subject (and not the simulated agent, beause we want to follow the exact learning path experienced by the subjects)
                        // after receiving feedback, update sv[t + 1]
                        sv[cue_gng[n,w,bl, t]] += ep_gng[n,w] * (rho_gng[n,w] * outcome_gng[n,w,bl,t] - sv[cue_gng[n,w,bl,t]]);

                        // update action values (based on the simulated response and the stochastic outcome we drew for this simulated trial)
                        if (pressed_gng[n,w,bl,t]) { // update go value
                            qv_g[cue_gng[n,w,bl, t]] += ep_gng[n,w] * (rho_gng[n,w] * outcome_gng[n,w,bl,t] - qv_g[cue_gng[n,w,bl,t]]);
                        }  
                        else { // update no-go value
                            qv_ng[cue_gng[n,w,bl,t]] += ep_gng[n,w] * (rho_gng[n,w] * outcome_gng[n,w,bl,t] - qv_ng[cue_gng[n,w,bl, t]]);
                        }
                    } // End of trials

                } // End of 3 blocks
            }                    
        }
    }
 }


