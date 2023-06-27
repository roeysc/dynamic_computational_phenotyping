// This model was adopted from the hBayesDM Package - https://github.com/CCS-Lab/hBayesDM


functions {

    real partial_sum(int[,,] choice_itc, int start, int end,  real[,,,] outcome_gng, int[,] idx_gng_obs,
                    int[,,,] pressed_gng, int[,,,] cue_gng,  int[,,] Tr_gng, int Bl, real[,] b_gng, real[,] pi_gng, 
                    real[,] xi_gng, real[,] ep_gng, real[,] rho_rew_pun_gng, real[,] rho_neut_gng, int W) {

	real lt = 0;
    real outcome;
	real rho;
	
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
                            // By default, take the true outcome that subjects received
                            outcome = outcome_gng[s,w,b,t];
                            rho = rho_rew_pun_gng[s,w];
                            // Only change the outcome and rho to use if this was a neutral outcome of 0
                            if ((outcome_gng[s,w,b,t]==0) && (cue_gng[s,w,b,t]<=2)){ //In Go2Win and NoGo2Win, neutral reward is a punishment
                                outcome = -1.0;
                                rho = rho_neut_gng[s,w];
                            } 
                            if ((outcome_gng[s,w,b,t]==0) && (cue_gng[s,w,b,t]>=3)){ //In Go2Avoid and NoGo2Avoid, neutral reward is a reward
                                outcome = 1.0;
                                rho = rho_neut_gng[s,w];
                            } 
                            
	                        sv[cue_gng[s,w,b,t]] += ep_gng[s,w] * (rho * outcome - sv[cue_gng[s,w,b,t]]);

	                        // update action values
	                        if (pressed_gng[s,w,b,t] == 1) { // update go value
	                            qv_g[cue_gng[s,w,b,t]] += ep_gng[s,w] * (rho * outcome - qv_g[cue_gng[s,w,b,t]]);
                            }
	                        else if  (pressed_gng[s,w,b,t] == 0) { // update no-go value
                                qv_ng[cue_gng[s,w,b,t]] += ep_gng[s,w] * (rho * outcome - qv_ng[cue_gng[s,w,b,t]]);
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

    // Exogenous variables
    int exo_q_num;                                                  // number of exogenous survey questions
    real U[N,W,exo_q_num];                                       // exogenous survey questions - missing weeks were linearly interpolated outside of Stan  

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
    int num_par = P_gng+1; // Because rho_gng is split into rho_rew_pun_gng and rho_neut
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];

    real xi_gng_pr[N,W];                 // gng, noise
    real ep_gng_pr[N,W];                 // gng, learning rate
    real rho_rew_pun_gng_pr[N,W];                // gng, rho
    real rho_neut_gng_pr[N,W];                // gng, rho
    real b_gng_pr[N,W];
    real pi_gng_pr[N,W];

}

transformed parameters {
    real b_gng[N,W];
    real pi_gng[N,W];
    real<lower=0, upper=1> xi_gng[N,W];                 // gng, noise
    real<lower=0, upper=1> ep_gng[N,W];                 // gng, learning rate
    real<lower=0> rho_rew_pun_gng[N,W];                // gng, rho
    real<lower=0> rho_neut_gng[N,W];                // gng, rho

    for (n in 1:N) {
        for (w in 1:W) {
            xi_gng[n,w] = Phi(mu_pr_sub[1,n] + sigma_pr_r[1]*xi_gng_pr[n,w]);
            ep_gng[n,w] = Phi(mu_pr_sub[2,n] + sigma_pr_r[2]*ep_gng_pr[n,w]);
            b_gng[n,w] = (mu_pr_sub[3,n] + sigma_pr_r[3]*b_gng_pr[n,w]);
            pi_gng[n,w] = exp(mu_pr_sub[4,n] + sigma_pr_r[4]*pi_gng_pr[n,w]);
            rho_rew_pun_gng[n,w] = exp(mu_pr_sub[5,n] + sigma_pr_r[5]*rho_rew_pun_gng_pr[n,w]);
            rho_neut_gng[n,w] = exp(mu_pr_sub[6,n] + sigma_pr_r[6]*rho_neut_gng_pr[n,w]);
        }
    }

}

model {
    // hyper parameters
    mu_pr[1]  ~ normal(0, 1.0); // We don't use the value from literature (xi_gng, irreducible noise)
    mu_pr[2]  ~ normal(0, 1.0); // We don't use the value from literature (ep_gng)
    mu_pr[3]  ~ normal(0, 1.0); // We don't use the value from literature  (b_gng, bias)
    mu_pr[4]  ~ normal(0, 1.0); // We don't use the value from literature  (pi_gng, Pavlovian bias)
    mu_pr[5]  ~ normal(0, 1.0); // We don't use the value from literature      (rho_gng, reward/punishmemt sensitivity)
    mu_pr[6]  ~ normal(0, 1.0); // We don't use the value from literature      (rho_gng, neutral outcome sensitivity)
    
    sigma_pr ~ normal(0, 1);

    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr[p], sigma_pr[p]);
    }

    to_vector(sigma_pr_r) ~ normal(0, 1);

    // individual parameters w/ Matt trick
    to_vector(to_matrix(xi_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(ep_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(rho_rew_pun_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(rho_neut_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(b_gng_pr)) ~ normal(0, 1);
    to_vector(to_matrix(pi_gng_pr)) ~ normal(0, 1);

    	target += reduce_sum(partial_sum, choice_itc, 1, outcome_gng, idx_gng_obs, pressed_gng, cue_gng, Tr_gng, Bl,
    	                    b_gng, pi_gng, xi_gng, ep_gng, rho_rew_pun_gng, rho_neut_gng, W);
}    



//////////////////////////////////////////////////////


 generated quantities {
    // For posterior predictive check
    real choice_ppc_gng[Tr_max_gng, Bl, W,N];
    real pGo_ppc_gng[Tr_max_gng, Bl, W,N];

    real log_lik_all_subs[N, W];

    real rho;
    real outcome;

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

	                for (b in 1:Bl) {

	                    vector[4] wv_g  = rep_vector(0.0, 4);  // action weight for go
	                    vector[4] wv_ng = rep_vector(0.0, 4); // action weight for nogo
	                    vector[4] qv_g = rep_vector(0.0, 4);  // Q value for go
	                    vector[4] qv_ng = rep_vector(0.0, 4); // Q value for nogo
	                    vector[4] sv = rep_vector(0.0, 4);    // stimulus value
	                    vector[4] pGo = rep_vector(0.0, 4);   // prob of go (press)

	                    for (t in 1:Tr_gng[n,w,b]) {
	                        wv_g[cue_gng[n,w,b,t]]  = qv_g[cue_gng[n,w,b,t]] + b_gng[n,w] + pi_gng[n,w] * sv[cue_gng[n,w,b,t]];
	                        wv_ng[cue_gng[n,w,b,t]] = qv_ng[cue_gng[n,w,b,t]];  // qv_ng is always equal to wv_ng (regardless of action)
	                        pGo[cue_gng[n,w, b,t]]   = inv_logit(wv_g[cue_gng[n,w, b,t]] - wv_ng[cue_gng[n,w,b,t]]);
                            pGo[cue_gng[n,w,b,t]]   *= (1 - xi_gng[n,w]);
                            pGo[cue_gng[n,w,b,t]]   += xi_gng[n,w]/2;
                            
                         // Keep track of the pGo value in this trial
                        pGo_ppc_gng[t,b,w,n] = pGo[cue_gng[n,w,b,t]];

                        // And keep track of the predicted choice made in this trial
                        choice_ppc_gng[t,b,w,n] = bernoulli_rng(pGo[cue_gng[n,w,b, t]]);

                        log_lik_all_subs[n,w] += bernoulli_lpmf(pressed_gng[n,w, b,t] | pGo[cue_gng[n,w,b, t]]);


	                        // after receiving feedback, update sv[t + 1]
                            // By default, take the true outcome that subjects received
                            outcome = outcome_gng[n,w,b,t];
                            rho = rho_rew_pun_gng[n,w];
                            // Only change the outcome and rho to use if this was a neutral outcome of 0
                            if ((outcome_gng[n,w,b,t]==0) && (cue_gng[n,w,b,t]<=2)){ //In Go2Win and NoGo2Win, neutral reward is a punishment
                                outcome = -1.0;
                                rho = rho_neut_gng[n,w];
                            } 
                            if ((outcome_gng[n,w,b,t]==0) && (cue_gng[n,w,b,t]>=3)){ //In Go2Avoid and NoGo2Avoid, neutral reward is a reward
                                outcome = 1.0;
                                rho = rho_neut_gng[n,w];
                            } 
                            
                            
                            sv[cue_gng[n,w,b,t]] += ep_gng[n,w] * (rho * outcome - sv[cue_gng[n,w,b,t]]);
	                        // update action values
	                        if (pressed_gng[n,w,b,t] == 1) { // update go value
                                qv_g[cue_gng[n,w,b,t]] += ep_gng[n,w] * (rho * outcome - qv_g[cue_gng[n,w,b,t]]);
                            }
	                        else if  (pressed_gng[n,w,b,t] == 0) { // update no-go value
                                qv_ng[cue_gng[n,w,b,t]] += ep_gng[n,w] * (rho * outcome - qv_ng[cue_gng[n,w,b,t]]);
                            }
	                    } // end of t loop
	                } // end of b loop
            }                    
        }
    }
 }
