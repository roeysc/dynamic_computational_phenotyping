

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
    int choice_itc[N,W,Tr_max_itc];     // choice itc - 0 for instant reward, 1 for delayed reward - subj x weeks x trials

}

transformed data {
    int num_par = P_itc;
    int exo_q_num_use = 2; // Number of exo effects to use
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];
    vector[num_par] pow_a_pr; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
    vector[num_par] mu_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters
    vector<lower=0>[num_par] sigma_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters



    real itc_k_pr[N,W];               // itc
    real itc_beta_pr[N,W];               // itc
    real pow_b_pr[N,num_par]; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions)

    matrix[num_par, exo_q_num_use] A_pr; //EXOMODEL
}


transformed parameters {
    real<lower=0> itc_k[N,W];        // itc
    real<lower=0> itc_beta[N,W];       // itc
    real<lower=0> pow_b[N,num_par]; // For modelling param = param_baseline  + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
    real pow_a_all_subjects[num_par,N];
    matrix[num_par, exo_q_num_use] A_all_subjects[N]; //EXOMODEL
        


    for (n in 1:N) {
        for (p in 1:num_par){
            pow_b[n,p] = exp(-2) + (exp(1)-exp(-2))*Phi(mu_pow_b_pr[p] + sigma_pow_b_pr[p]*pow_b_pr[n,p]); // Constrain b to values that don't lead to flat learning curves (which effectively remove the learning from the learning term, causing identifiability issues)
        }
        real ses_count = 1.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
        int first_session; // Boolean
        first_session = 1;
        
        // Create the subject-specific learning term pow_a
        vector[num_par] pow_a;
        pow_a[1] = (2.0*Phi_approx(pow_a_pr[1])-1.0).*mu_pr_sub[1,n]; // each pow_a_frac is a fraction of the corresponding mu_pr_sub (multiplied by 2 to allow for greater learning)
        pow_a[2] = (2.0*Phi_approx(pow_a_pr[2])-1.0).*mu_pr_sub[2,n]; // each pow_a_frac is a fraction of the corresponding mu_pr_sub (multiplied by 2 to allow for greater learning)
        pow_a_all_subjects[1,n] = pow_a[1];
        pow_a_all_subjects[2,n] = pow_a[2];

        // Create the subject-specific exogenous term
        matrix[num_par, exo_q_num_use] A; //EXOMODEL
        for (p in 1:num_par){
            for (q in 1:exo_q_num_use){
                A[p,q] = (2.0*Phi_approx(A_pr[p,q])-1.0).*mu_pr_sub[p,n]; // each exogenous loading is a fraction of the corresponding mu_pr
            }
        }
        A_all_subjects[n] = A;

        //real prev_itc_k; // Accumulated noise in the unconstrained space
        //real prev_itc_beta; // Accumulated noise in the unconstrained space
        
        for (w in 1:W) {
            itc_k[n,w] = 0; // Initialize to avoid parameters with undefined values
            itc_beta[n,w] = 0;
            if (idx_itc_obs[n,w] != 0) { // if week exists
                //if (first_session == 1) {
                //    prev_itc_k = 0;
                //    prev_itc_beta = 0;
                //    first_session = 0;
                //}
                itc_k[n,w] = exp((mu_pr_sub[1,n])   +   (A[1,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1])) + sigma_pr_r[1]*itc_k_pr[n,w]);
                itc_beta[n,w] = exp((mu_pr_sub[2,n])   +(A[2,] * to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2])) + sigma_pr_r[2]*itc_beta_pr[n,w]);
            
                //prev_itc_k =    prev_itc_k + sigma_pr_r[1]*itc_k_pr[n,w]; // Accumulate noise over time to model a slow diffusing process
                //prev_itc_beta = prev_itc_beta + sigma_pr_r[2]*itc_beta_pr[n,w];
                ses_count = ses_count+1.0;
            }
        }
    }

}

model {
    // hyper parameters
    mu_pr[1]  ~ normal(-4.17, 0.65); // literature based prior for the k parameter
    mu_pr[2]  ~ normal(-1.4,1); // literature based prior for the beta parameter

    sigma_pr ~ normal(0, 1.0); // We could shrink this, because we assume subject are overall similar to each other
    
    pow_a_pr ~ normal(0, 1.0); // We could shrink this, because we want dynamical parameters to be close to zero
    mu_pow_b_pr  ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. We do shrink pow_a (C in the overlead notation), so I think that's enough
  

    to_vector(A_pr) ~ normal(0, 1.0); 
    
    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr[p], sigma_pr[p]);
    }

    to_vector(sigma_pr_r) ~ normal(0, 1);
    sigma_pow_b_pr ~ normal(0, 1);

    // individual parameters w/ Matt trick
    to_vector(to_matrix(itc_k_pr)) ~ normal(0, 1);
    to_vector(to_matrix(itc_beta_pr)) ~ normal(0, 1);
    // same for learning curves
    to_vector(to_matrix(pow_b_pr)) ~ normal(0, 1);

    	target += reduce_sum(partial_sum, choice_itc, 1, delay_later, amount_later, amount_sooner, itc_k, itc_beta, 
    	                    idx_itc_obs, Tr_itc, W);
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
        int first_session; // Boolean
        first_session = 1;
        
        // Get the subject-specific learning term pow_a
        vector[num_par] pow_a;
        pow_a = to_vector(pow_a_all_subjects[,n]);

        // Get the subject-specific exogenous term
        matrix[num_par, exo_q_num_use] A; //EXOMODEL
        A = A_all_subjects[n];
        
        real prev_itc_k; // Previous value in unconstrained space
        real prev_itc_beta; // Previous value in unconstrained space

        for (w in 1:W) {
            if (idx_itc_obs[n,w] != 0) { // if week exists
                //if (first_session == 1) {
                //    prev_itc_k = 0;
                //    prev_itc_beta = 0;
                //    first_session = 0;
                //}
                
                
                // By using element-wise multiplication we make sure to have both A*U and (-A)*U, each with its own weight
                //prev_itc_k =     w_mu[1]*(mu_pr_sub[1,n])   +   w_prv[1]*(prev_itc_k)   +   w_exo[1]*(mu_pr_sub[1,n] + A[1,] * to_vector(U[n,w,]))   +    w_lrn[1]*(mu_pr_sub[1,n] - pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1])) + sigma_pr_r[1]*itc_k_pr[n,w];
                //prev_itc_beta =     w_mu[2]*(mu_pr_sub[2,n])   +   w_prv[2]*(prev_itc_beta)   +   w_exo[2]*(mu_pr_sub[2,n] + A[2,] * to_vector(U[n,w,]))   +    w_lrn[2]*(mu_pr_sub[2,n] - pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2])) + sigma_pr_r[2]*itc_beta_pr[n,w];
                
                // Fill in the terms we will use to calculate the variance explained of each term                
                prev[1,n,w] = 0;
                exo1[1,n,w] = (A[1,1]*(U[n,w,1]) );
                exo2[1,n,w] = (A[1,2]*(U[n,w,2]) );
                lrn[1,n,w] = (-pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1]));
                noise[1,n,w] = sigma_pr_r[1]*itc_k_pr[n,w];

                prev[2,n,w] = 0;
                exo1[2,n,w] = (A[2,1]*(U[n,w,1]) );
                exo2[2,n,w] = (A[2,2]*(U[n,w,2]) );
                lrn[2,n,w] = (-pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2]));
                noise[2,n,w] = sigma_pr_r[2]*itc_beta_pr[n,w];

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
            A_frac[p,q] = (2.0*Phi_approx(A_pr[p,q])-1.0)  * sgn;
        }
    }

    
}
 