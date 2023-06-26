// In this version, the evidence accumulation rate delta is modeled with an exp() transform. This ensures positivity, and also allows obtaining high values more easily (indeed, we expect the raw delta to be quite high, on the order of 50-100, since we then multiply it by the coherence fraction to obtain the trial-specific evidence accumulation rate).

functions {

    real partial_sum(int[,,] choice_itc, int start, int end, real[,,] RTu_rdm, real[,,] RTl_rdm, real[,,] Cohu_rdm, 
                    real[,,] Cohl_rdm, real[,] delta_rdm, real[,] alpha_rdm, real[,] tau_rdm, int[,] idx_rdm_obs, 
                    int[,] Nu_rdm, int[,] Nl_rdm, int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {        	    			

 	        if (idx_rdm_obs[s,w] != 0) { 	 
	                vector[Nu_rdm[s,w]] delta_cohu = delta_rdm[s,w]*to_vector(Cohu_rdm[s,w,:Nu_rdm[s,w]]);
	                vector[Nl_rdm[s,w]] delta_cohl = delta_rdm[s,w]*to_vector(Cohl_rdm[s,w,:Nl_rdm[s,w]]);

	                lt += wiener_lpdf(RTu_rdm[s,w,:Nu_rdm[s,w]] | alpha_rdm[s,w], tau_rdm[s,w], 0.5, delta_cohu);
                    lt += wiener_lpdf(RTl_rdm[s,w,:Nl_rdm[s,w]] | alpha_rdm[s,w], tau_rdm[s,w], 0.5, -delta_cohl); 
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

    // Random dot motion
    int P_rdm;
    int<lower=0> idx_rdm_obs[N,W];                            // Indices of weeks WITH data
    int<lower=0> Nu_max_rdm;                            // Max (across subjects) number of upper boundary responses
    int<lower=0> Nl_max_rdm;                        // Max (across subjects) number of lower boundary responses
    int<lower=0> Nu_rdm[N,W];                     // Number of upper boundary responses for each subj
    int<lower=0> Nl_rdm[N,W];                         // Number of lower boundary responses for each subj
    real RTu_rdm[N, W, Nu_max_rdm];                // upper boundary response times
    real RTl_rdm[N, W, Nl_max_rdm];                    // lower boundary response times
    real Cohu_rdm[N, W, Nu_max_rdm];                           // coherence for correct trials
    real Cohl_rdm[N, W, Nl_max_rdm];                   // coherence for incorrect trials
    matrix[N,W] minRT_rdm;                          // minimum RT for each subject of the observed data
    real RTbound_rdm;                        // choice nc - subj x weeks x trials
}

transformed data {
    int num_par = P_rdm;
    int num_par_lrn = num_par-1; // Number of parameters for which we model learning
    int num_par_exo = num_par-1; // Number of parameters for which we model learning
    int exo_q_num_use = 2; // How many exo loadings to actually use
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];
    vector[num_par_lrn] pow_a_pr; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
    vector[num_par_lrn] mu_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters
    vector<lower=0>[num_par_lrn] sigma_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters

    //vector[num_par] w_prv_pr; // Weight for the previous week in the slow-dynamics model (one parameter shared across all subjects)
    

    real pow_b_pr[N,num_par_lrn]; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions)

    matrix[num_par_exo, exo_q_num_use] A_pr; //EXOMODEL

    real alpha_rdm_pr[N,W];              // rdm
    real delta_rdm_pr[N,W];              // rdm
    real tau_rdm_pr[N,W];                // rdm

}

transformed parameters {
    real tau_rdm[N,W];
    real<lower=0> alpha_rdm[N,W];
    real<lower=0> delta_rdm[N,W];
    real pow_a_all_subjects[num_par_lrn,N];
    matrix[num_par_exo, exo_q_num_use] A_all_subjects[N]; //EXOMODEL
    real<lower=0> pow_b[N,num_par_lrn]; // For modelling param = param_baseline  + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
  
    //vector<lower=0>[num_par] w_prv; // Weight for the previous week in the slow-dynamics model (one parameter shared across all subjects)
    
    for (n in 1:N) {
        for (p in 1:num_par_lrn){
            pow_b[n,p] = exp(-2) + (exp(1)-exp(-2))*Phi(mu_pow_b_pr[p] + sigma_pow_b_pr[p]*pow_b_pr[n,p]); // Constrain b to values that don't lead to flat learning curves (which effectively remove the learning from the learning term, causing identifiability issues)
        }
        real ses_count = 1.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
        int first_session; // Boolean
        first_session = 1;
        
         // Create the subject-specific learning term pow_a
        vector[num_par_lrn] pow_a;
        for (p in 1:2){ // ROEY: IF THIS FAILS, SEPARATE THE LOOP TO 1 AND 2, AND FOR PARAM2 USE 1.0*PHI_APPROX INSTEAD OF 2.0*PHI_APPROX
            pow_a[p] = (2.0*Phi_approx(pow_a_pr[p])-1.0).*mu_pr_sub[p,n]; // each pow_a_frac is a fraction of the corresponding mu_pr_sub (multiplied by 2 to allow for greater learning)
            pow_a_all_subjects[p,n] = pow_a[p];
        }
        pow_a = to_vector(pow_a_all_subjects[,n]);
        
        // Create the subject-specific exogenous term
        matrix[num_par_exo, exo_q_num_use] A; //EXOMODEL
        for (p in 1:num_par_exo){
            for (q in 1:exo_q_num_use){
                A[p,q] = (2.0*Phi_approx(A_pr[p,q])-1.0).*mu_pr_sub[p,n]; // each exogenous loading is a fraction of the corresponding mu_pr
            }
        }
        A_all_subjects[n] = A;
        
        

        //real prev_alpha_rdm; // Previous value in unconstrained space
        //real prev_delta_rdm; // Previous value in unconstrained space
        //real prev_tau_rdm; // Previous value in unconstrained space
        
        for (w in 1:W) {
            alpha_rdm[n,w] = 0; // Initialize to avoid parameters with undefined values
            delta_rdm[n,w] = 0;
            tau_rdm[n,w] = 0; // Initialize to avoid parameters with undefined values
            if (idx_rdm_obs[n,w] != 0) { // if week exists
                //if (first_session == 1) {
                //    prev_alpha_rdm = mu_pr_sub[1,n];
                //    prev_delta_rdm = mu_pr_sub[2,n];
                //    prev_tau_rdm = mu_pr_sub[3,n];
                //    first_session = 0;
                //}
                alpha_rdm[n,w] =        exp((mu_pr_sub[1,n])   +  (A[1,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (-  pow_a[1] + pow_a[1] *   pow(ses_count,-pow_b[n,1])) + sigma_pr_r[1]*alpha_rdm_pr[n,w]);
                delta_rdm[n,w] =        exp((mu_pr_sub[2,n])   +  (A[2,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (-  pow_a[2] + pow_a[2] *   pow(ses_count,-pow_b[n,2])) + sigma_pr_r[2]*delta_rdm_pr[n,w]);
                tau_rdm[n,w] =    inv_logit((mu_pr_sub[3,n])   + sigma_pr_r[3]*tau_rdm_pr[n,w]) * (minRT_rdm[n,w] - RTbound_rdm - 0.0001) + RTbound_rdm;   // TauBound

                //prev_alpha_rdm =  (w_mu[1]*(mu_pr_sub[1,n])   +   w_prv[1]*(prev_alpha_rdm)   +  w_exo[1]*(mu_pr_sub[1,n] +   A[1,] * to_vector(U[n,w,]))   +    w_lrn[1]*(mu_pr_sub[1,n] -  pow_a[1] + pow_a[1] *   pow(ses_count,-pow_b[n,1])) + sigma_pr_r[1]*alpha_rdm_pr[n,w]);
                //prev_delta_rdm =  (w_mu[2]*(mu_pr_sub[2,n])   +   w_prv[2]*(prev_delta_rdm)   +  w_exo[2]*(mu_pr_sub[2,n] +   A[2,] * to_vector(U[n,w,]))   +    w_lrn[2]*(mu_pr_sub[2,n] -  pow_a[2] + pow_a[2] *   pow(ses_count,-pow_b[n,2])) + sigma_pr_r[2]*delta_rdm_pr[n,w]);
                //prev_tau_rdm =    (w_mu[3]*(mu_pr_sub[3,n])   +   w_prv[3]*(prev_tau_rdm)   +    w_exo[3]*(mu_pr_sub[3,n] +   A[3,] * to_vector(U[n,w,]))   +    w_lrn[3]*(mu_pr_sub[3,n] -  pow_a[3] + pow_a[3] *   pow(ses_count,-pow_b[n,3])) + sigma_pr_r[3]*tau_rdm_pr[n,w]);
                
                ses_count = ses_count+1.0;
         
            }
        }
    }

}

model {
    // hyper parameters
    mu_pr[1]  ~ normal(0, 1.0); // We use the value from literature
    mu_pr[2]  ~ normal(0, 1.0); // We use the value from literature
    mu_pr[3]  ~ normal(0, 1.0); // We use the value from literature, but increased the STD by x100, because the choice of minRT to use was quite arbitrary (we used the 5th percentile across the population, but in practice we use the minRT of each session separtely)
    sigma_pr ~ normal(0, 1.0); // We could shrink this, because we assume subject are overall similar to each other
    
    pow_a_pr ~ normal(0, 1.0); // We could shrink this, because we want dynamical parameters to be close to zero
    mu_pow_b_pr  ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. We do shrink pow_a (C in the overlead notation), so I think that's enough
    to_vector(to_matrix(pow_b_pr)) ~ normal(0, 1.0);

    //w_prv_pr ~ normal(0, 1.0);
    
    to_vector(A_pr) ~ normal(0, 1.0); // ***  We shrink this, because we want dynamical parameters to be close to zero
   
    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr[p], sigma_pr[p]);
    }

    to_vector(sigma_pr_r) ~ normal(0, 1.0); // Here we DON'T shrink the week-to-week noise
    sigma_pow_b_pr ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. &&&TODO: MAYBE WE DO WANT TO SHRINK THIS UNDER THE ASSUMPTION THAT SUBJECTS ARE SIMILAR TO EACH OTHER?

    to_vector(to_matrix(alpha_rdm_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(delta_rdm_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(tau_rdm_pr)) ~ normal(0, 1.0);

    	target += reduce_sum(partial_sum, choice_itc, 1, RTu_rdm, RTl_rdm, Cohu_rdm, Cohl_rdm, delta_rdm,alpha_rdm, 
    	                    tau_rdm, idx_rdm_obs, Nu_rdm, Nl_rdm, W);
}  



generated quantities {
    // Keep record of each of the dynamical terms, so we can calculate their variance explained
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
             }
        }
   }


   for (n in 1:N) {
        real ses_count = 1.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
        
        // Get the subject-specific learning term pow_a
        vector[num_par-1] pow_a;
        pow_a = to_vector(pow_a_all_subjects[,n]);

        // Get the subject-specific exogenous term
        matrix[num_par-1, exo_q_num_use] A; //EXOMODEL
        A = A_all_subjects[n];
        
        for (w in 1:W) {
            if (idx_rdm_obs[n,w] != 0) { // if week exists
                // Fill in the terms we will use to calculate the variance explained of each term                
                exo1[1,n,w] = A[1,1]*(U[n,w,1]);
                exo2[1,n,w] = A[1,2]*(U[n,w,2]);
                lrn[1,n,w] = -pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1]);
                noise[1,n,w] = sigma_pr_r[1]*alpha_rdm_pr[n,w];

                exo1[2,n,w] = A[2,1]*(U[n,w,1]);
                exo2[2,n,w] = A[2,2]*(U[n,w,2]);
                lrn[2,n,w] = -pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2]);
                noise[2,n,w] = sigma_pr_r[2]*delta_rdm_pr[n,w];

                ses_count = ses_count+1.0;
            }
        
        }
    }
    
    // *** Get the general exo and learning terms, to calculate their probability of direction
    vector[num_par-1] pow_a_frac;
    matrix[num_par-1,exo_q_num_use] A_frac;
    real sgn; // we multiply the terms by the sign of the parameter's baseline so we know of the effect is increasing or decreasing the parameter (for leaerning, s negative sign means increasing - checl out the opwerpower law formulation to see why) 
    for (p in 1:num_par-1){
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