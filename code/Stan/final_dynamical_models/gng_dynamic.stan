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
    int<lower=-1, upper=1> choice_itc[N,W,Tr_max_itc];     // choice itc - 0 for instant reward, 1 for delayed reward - subj x weeks x trials

    int P_gng;
    int<lower=0> idx_gng_obs[N,W];                            // Indices of weeks WITH data
    int<lower=0> Bl;
    int<lower=1> Tr_max_gng;
    int<lower=0, upper=Tr_max_gng> Tr_gng[N,W,Bl];
    int<lower=0, upper=4> cue_gng[N,W,Bl,Tr_max_gng];
    int pressed_gng[N,W,Bl, Tr_max_gng];
    real outcome_gng[N, W, Bl, Tr_max_gng];  

}

transformed data {
    int num_par = P_gng + 1; // Because rho_gng is split into rho_rew_pun_gng and rho_neut
    int num_par_lrn = P_gng + 1; // Number of parameters for which we model learning.
    int exo_q_num_use = 2; // Number of exo effects to use
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];
    vector[num_par_lrn] pow_a_pr; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
    vector[num_par_lrn] mu_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters
    vector<lower=0>[num_par_lrn] sigma_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters
    
    real pow_b_pr[N,num_par_lrn]; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions)

    matrix[num_par, exo_q_num_use] A_pr; //EXOMODEL
  
    real xi_gng_pr[N,W];                 // gng, noise
    real ep_gng_pr[N,W];                 // gng, learning rate
    real rho_rew_pun_gng_pr[N,W];                // gng, rho
    real rho_neut_gng_pr[N,W];                // gng, rho
    real b_gng_pr[N,W];
    real pi_gng_pr[N,W];

}

transformed parameters {
    real<lower=0, upper=1> xi_gng[N,W];                 // gng, noise
    real<lower=0, upper=1> ep_gng[N,W];                 // gng, learning rate
    real b_gng[N,W];
    real pi_gng[N,W];
    real<lower=0> rho_rew_pun_gng[N,W];                // gng, rho
    real<lower=0> rho_neut_gng[N,W];                // gng, rho
    real pow_a_all_subjects[num_par_lrn,N];
    matrix[num_par, exo_q_num_use] A_all_subjects[N]; //EXOMODEL
    real<lower=0> pow_b[N,num_par_lrn]; // For modelling param = param_baseline  + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training


    for (n in 1:N) {
        for (p in 1:num_par_lrn){
            pow_b[n,p] = exp(-2) + (exp(1)-exp(-2))*Phi(mu_pow_b_pr[p] + sigma_pow_b_pr[p]*pow_b_pr[n,p]); // Constrain b to values that don't lead to flat learning curves (which effectively remove the learning from the learning term, causing identifiability issues)
        }
        real ses_count = 1.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
        int first_session; // Boolean
        first_session = 1;
        
         // Create the subject-specific learning term pow_a
        vector[num_par_lrn] pow_a_frac;
        vector[num_par_lrn] pow_a;
        for (p in 1:num_par_lrn){
            pow_a_frac[p] = (2.0*Phi_approx(pow_a_pr[p])-1.0).*mu_pr_sub[p,n]; // each pow_a_frac is a fraction of the corresponding mu_pr_sub (multiplied by 2 to allow for greater learning)
            pow_a_all_subjects[p,n] = pow_a_frac[p];
        }
        pow_a = to_vector(pow_a_all_subjects[,n]);
        
        // Create the subject-specific exogenous term
        matrix[num_par, exo_q_num_use] A; //EXOMODEL
        for (p in 1:num_par){
            for (q in 1:exo_q_num_use){
                A[p,q] = (2.0*Phi(A_pr[p,q])-1.0).*mu_pr_sub[p,n]; // each exogenous loading is a fraction of the corresponding mu_pr
            }
        }
        A_all_subjects[n] = A;
        
        for (w in 1:W) {
            xi_gng[n,w] = 0; // Initialize to avoid parameters with undefined values
            ep_gng[n,w] = 0;
            b_gng[n,w] = 0; // Initialize to avoid parameters with undefined values
            pi_gng[n,w] = 0;
            rho_rew_pun_gng[n,w] = 0; // Initialize to avoid parameters with undefined values
            rho_neut_gng[n,w] = 0;
            if (idx_gng_obs[n,w] != 0) { // if week exists
                xi_gng[n,w] =  Phi((mu_pr_sub[1,n])   +  (A[1,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1])) + sigma_pr_r[1]*xi_gng_pr[n,w]);
                ep_gng[n,w] =  Phi((mu_pr_sub[2,n])   +  (A[2,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2])) + sigma_pr_r[2]*ep_gng_pr[n,w]);
                b_gng[n,w] =     ((mu_pr_sub[3,n])    +  (A[3,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[3] + pow_a[3] * pow(ses_count,-pow_b[n,3])) + sigma_pr_r[3]*b_gng_pr[n,w]);
                pi_gng[n,w] =    exp((mu_pr_sub[4,n])    +  (A[4,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[4] + pow_a[4] * pow(ses_count,-pow_b[n,4])) + sigma_pr_r[4]*pi_gng_pr[n,w]);
                rho_rew_pun_gng[n,w] = exp((mu_pr_sub[5,n])   +  (A[5,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[5] + pow_a[5] * pow(ses_count,-pow_b[n,5])) + sigma_pr_r[5]*rho_rew_pun_gng_pr[n,w]);
                rho_neut_gng[n,w] = exp((mu_pr_sub[6,n])   +  (A[6,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[6] + pow_a[6] * pow(ses_count,-pow_b[n,6])) + sigma_pr_r[6]*rho_neut_gng_pr[n,w]);

                ses_count = ses_count+1.0;
            }
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

    sigma_pr ~ normal(0, 1.0); // We could shrink this, because we assume subject are overall similar to each other

    
    pow_a_pr ~ normal(0, 1.0); // We could shrink this, because we want dynamical parameters to be close to zero
    mu_pow_b_pr  ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. We do shrink pow_a (C in the overlead notation), so I think that's enough
    to_vector(to_matrix(pow_b_pr)) ~ normal(0, 1.0);

    to_vector(A_pr) ~ normal(0, 1); // ***  We shrink this, because we want dynamical parameters to be close to zero
    
    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr[p], sigma_pr[p]);
    }

    to_vector(sigma_pr_r) ~ normal(0, 1.0); // Here we DON'T shrink the week-to-week noise
    sigma_pow_b_pr ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. &&&TODO: MAYBE WE DO WANT TO SHRINK THIS UNDER THE ASSUMPTION THAT SUBJECTS ARE SIMILAR TO EACH OTHER?

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






generated quantities {
    // Keep record of each of the dynamical terms, so we can calculate their variance explained
    real prev[num_par,N,W];
    real exo1[num_par,N,W];
    real exo2[num_par,N,W];
    real lrn[num_par,N,W];
    real noise[num_par,N,W];
    
    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {
            for (p in 1:num_par) {
                 lrn[p,n,w]=-2;
                 exo1[p,n,w]=-2;
                 exo2[p,n,w]=-2;
                 noise[p,n,w]=-2;
                 prev[p,n,w]=-2;
             }
        }
   }


   for (n in 1:N) {
        real ses_count = 1.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
        
        // Get the subject-specific learning term pow_a
        vector[num_par] pow_a;
        pow_a = to_vector(pow_a_all_subjects[,n]);

        // Get the subject-specific exogenous term
        matrix[num_par, exo_q_num_use] A; //EXOMODEL
        A = A_all_subjects[n];
        
        for (w in 1:W) {
            if (idx_gng_obs[n,w] != 0) { // if week exists

                // Fill in the terms we will use to calculate the variance explained of each term                
                for (p in 1:num_par) {
                    exo1[p,n,w] = A[p,1]*(U[n,w,1]);
                    exo2[p,n,w] = A[p,2]*(U[n,w,2]);
                    lrn[p,n,w] = (-pow_a[p] + pow_a[p] * pow(ses_count,-pow_b[n,p]));
                }
                noise[1,n,w] = sigma_pr_r[1]*xi_gng_pr[n,w];
                noise[2,n,w] = sigma_pr_r[2]*ep_gng_pr[n,w];
                noise[3,n,w] = sigma_pr_r[3]*b_gng_pr[n,w];
                noise[4,n,w] = sigma_pr_r[4]*pi_gng_pr[n,w];
                noise[5,n,w] = sigma_pr_r[5]*rho_rew_pun_gng_pr[n,w];
                noise[6,n,w] = sigma_pr_r[6]*rho_neut_gng_pr[n,w];
                
                ses_count = ses_count+1.0;
            }
        
        }
    }
    
    // *** Get the general exo and learning terms, to calculate their probability of direction
    vector[num_par] pow_a_frac;
    matrix[num_par,exo_q_num_use] A_frac;
    real sgn; // we multiply the terms by the sign of the parameter's baseline so we know of the effect is increasing or decreasing the parameter (for leaerning, s negative sign means increasing - checl out the opwerpower law formulation to see why) 
    for (p in 1:num_par){
       if (mu_pr[p]<0){
            sgn=-1;
        }else{
            sgn=1;
        }
        pow_a_frac[p] = (2.0*Phi_approx(pow_a_pr[p])-1.0) * sgn;
        for (q in 1:exo_q_num_use){
            A_frac[p,q] = (2.0*Phi_approx(A_pr[p,q])-1.0) * sgn;
        }
    }

}
