

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

    // Two-armed bandit
    int<lower=1> P_smb; // number of predictors
    int<lower=0> idx_smb_obs[N,W];                            // Indices of weeks WITH data
    int<lower=1> Tr_max_smb; // max number of trials
    int<lower=0, upper=Tr_max_smb> Tr_smb[N, W]; 
    real x_smb[N, W, Tr_max_smb, P_smb];  // predictor matrix
    int y_smb[N, W, Tr_max_smb];  // Outcome - sub X trials

    // Numerocity comparison 
    int P_nc;
    int<lower=1> W_nc_obs[N];                                          // Number of weeks WITH data
    int<lower=0> idx_nc_obs[N,W];                              // Indices of weeks WITH data
    int<lower=0> Tr_max_nc;                                          // Max number of trials across subjects across weeks
    int<lower=0> Tr_nc[N, W];                                // Number of trials for each subj for each week
    int deltaM[N, W, Tr_max_nc];                              // Mu of dots - subj x weeks x trials
    real TotalS[N, W, Tr_max_nc];                             // Sd of dots - subj x weeks x trials
    int choice_nc[N, W, Tr_max_nc];                           // choice nc - subj x weeks x trials

    int P_gng;
    int<lower=0> idx_gng_obs[N,W];                            // Indices of weeks WITH data
    int<lower=0> Bl;
    int<lower=1> Tr_max_gng;
    int<lower=0, upper=Tr_max_gng> Tr_gng[N,W,Bl];
    int<lower=0, upper=4> cue_gng[N,W,Bl,Tr_max_gng];
    int<lower=-1, upper=1> pressed_gng[N,W,Bl, Tr_max_gng];
    real outcome_gng[N, W, Bl, Tr_max_gng];  

    int P_lt;
    int<lower=0> idx_lt_obs[N,W];                              // Indices of weeks WITH data
    int<lower=0> Tr_max_lt;                                          // Max number of trials across subjects across weeks
    int<lower=0> Tr_lt[N, W];                                // Number of trials for each subj for each week
    real<lower=0> hi_p_lt[N, W,Tr_max_lt];
    real<lower=0> hi_narr_lt[N, W, Tr_max_lt];
    real<lower=0> lo_narr_lt[N, W, Tr_max_lt];
    real<lower=0> hi_wide_lt[N, W, Tr_max_lt];
    real<lower=0> lo_wide_lt[N, W, Tr_max_lt];
    int<lower=-1, upper=1> choice_lt[N, W, Tr_max_lt]; // Risky=0; Safe=1. Roey: Why is -1 an option here?   

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
    int num_par = P_smb;
    int exo_q_num_use = 2; // number of exo effects to use
}

parameters {
    vector[num_par] mu_pr;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];
    vector[num_par] pow_a_pr; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
    vector[num_par] mu_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters
    vector<lower=0>[num_par] sigma_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters
    
    real params_smb_pr[N, W, 3];            // smb 3 params

    real pow_b_pr[N,num_par]; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions)

    matrix[num_par, exo_q_num_use] A_pr; //EXOMODEL
}

transformed parameters {
    real params_smb[N, W, 3];           // smb
    real pow_a_all_subjects[num_par,N];
    matrix[num_par, exo_q_num_use] A_all_subjects[N]; //EXOMODEL
    real<lower=0> pow_b[N,num_par]; // For modelling param = param_baseline  + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
    

    for (n in 1:N) {
        for (p in 1:num_par){
            pow_b[n,p] = exp(-2) + (exp(1)-exp(-2))*Phi(mu_pow_b_pr[p] + sigma_pow_b_pr[p]*pow_b_pr[n,p]); // Constrain b to values that don't lead to flat learning curves (which effectively remove the learning from the learning term, causeing identifiability issues)
        }
        real ses_count = 1.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
        int first_session; // Boolean
        first_session = 1;
     
         // Create the subject-specific learning term pow_a
        vector[num_par] pow_a;
        for (p in 1:num_par){
            pow_a[p] = (2.0*Phi_approx(pow_a_pr[p])-1.0).*mu_pr_sub[p,n]; // each pow_a_frac is a fraction of the corresponding mu_pr_sub (multiplied by 2 to allow for greater learning)
            pow_a_all_subjects[p,n] = pow_a[p];
        }
        
        // Create the subject-specific exogenous term
        matrix[num_par, exo_q_num_use] A; //EXOMODEL
        for (p in 1:num_par){
            for (q in 1:exo_q_num_use){
                A[p,q] = (2.0*Phi_approx(A_pr[p,q])-1.0).*mu_pr_sub[p,n]; // each exogenous loading is a fraction of the corresponding mu_pr
            }
        }
        A_all_subjects[n] = A;
        
     
        //real prev_smb_wV; // Previous value in unconstrained space
        //real prev_smb_wsignV0verTU; // Previous value in unconstrained space
        //real prev_smb_wRU; // Previous value in unconstrained space


        for (w in 1:W) {
            params_smb[n,w,1] = 0; // Initialize to avoid parameters with undefined values
            params_smb[n,w,2] = 0;
            params_smb[n,w,3] = 0;
            if (idx_smb_obs[n,w] != 0) { // if week exists
                //if (first_session == 1) {
                //    prev_smb_wV = mu_pr_sub[1,n];
                //    prev_smb_wsignV0verTU = mu_pr_sub[2,n];
                //    prev_smb_wRU = mu_pr_sub[3,n];
                //    first_session = 0;
                //}
                params_smb[n,w,1] = (mu_pr_sub[1,n] )   +   (A[1,]*to_vector(U[n,w,1:exo_q_num_use]))   +  (- pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1])) + sigma_pr_r[1]*params_smb_pr[n,w,1];     // smb                                  
                params_smb[n,w,2] = (mu_pr_sub[2,n] )   +   (A[2,]*to_vector(U[n,w,1:exo_q_num_use]))   +  (- pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2])) + sigma_pr_r[2]*params_smb_pr[n,w,2];     // smb                                  
                params_smb[n,w,3] = (mu_pr_sub[3,n] )   +   (A[3,]*to_vector(U[n,w,1:exo_q_num_use]))   +  (- pow_a[3] + pow_a[3] * pow(ses_count,-pow_b[n,3])) + sigma_pr_r[3]*params_smb_pr[n,w,3];     // smb                                  

                //prev_smb_wV           = prev_smb_wV           + sigma_pr_r[1]*params_smb_pr[n,w,1];
                //prev_smb_wsignV0verTU = prev_smb_wsignV0verTU + sigma_pr_r[2]*params_smb_pr[n,w,2];
                //prev_smb_wRU          = prev_smb_wRU          + sigma_pr_r[3]*params_smb_pr[n,w,3];
                ses_count = ses_count+1.0;
                
            }
        }

    }

}

model {
    // hyper parameters
    // Phenotype priors are based on Gershman (2019) "Uncertainty and exploration" (Figure 3). In fact, we had to download the data from GitHub, add sign(V)/TU as a predictor, and set standardize=True to calculate the fixed effects estimates we used as priors here. The STD were calculated as SE*sqrt(46), as there were 46 subjects in that 2019 paper. 
    mu_pr[1]  ~ normal(0, 1.0); // Based on Sam's paper (but with sign(V)/TU - see https://docs.google.com/document/d/1lB45Ixnj40xqTFhZz17mXz3vAaX8Pv3y1CVZx6HqlXg/edit?usp=sharing)
    mu_pr[2]  ~ normal(0, 1.0); // Based on Sam's paper (but with sign(V)/TU - see https://docs.google.com/document/d/1lB45Ixnj40xqTFhZz17mXz3vAaX8Pv3y1CVZx6HqlXg/edit?usp=sharing)
    mu_pr[3]  ~ normal(0, 1.0); // Based on Sam's paper (but with sign(V)/TU - see https://docs.google.com/document/d/1lB45Ixnj40xqTFhZz17mXz3vAaX8Pv3y1CVZx6HqlXg/edit?usp=sharing)
    sigma_pr ~ normal(0, sqrt(1.0)); // We could shrink this, because we assume subject are overall similar to each other
    pow_a_pr ~ normal(0, sqrt(1.0)); // *** We could shrink this, because we want dynamical parameters to be close to zero
    mu_pow_b_pr  ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. We do shrink pow_a (C in the overlead notation), so I think that's enough
    to_vector(A_pr) ~ normal(0, sqrt(1.0)); // ***  We shrink this, because we want dynamical parameters to be close to zero
    
    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr[p], sigma_pr[p]);
    }

    to_vector(sigma_pr_r) ~ normal(0, 1.0);
    sigma_pow_b_pr ~ normal(0, 1.0);

    // individual parameters w/ Matt trick
    to_vector(to_matrix(params_smb_pr[,,1])) ~ normal(0,1.0);
    to_vector(to_matrix(params_smb_pr[,,2])) ~ normal(0,1.0);
    to_vector(to_matrix(params_smb_pr[,,3])) ~ normal(0,1.0);
    // same for learning curves
    to_vector(to_matrix(pow_b_pr)) ~ normal(0, 1);

    	target += reduce_sum(partial_sum, choice_itc, 1, y_smb, params_smb, idx_smb_obs, x_smb, Tr_smb, W);
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
            if (idx_smb_obs[n,w] != 0) { // if week exists

                // Fill in the terms we will use to calculate the variance explained of each term                
                for (p in 1:num_par) {
                    exo1[p,n,w] = (A[p,1]*(U[n,w,1]) );
                    exo2[p,n,w] = (A[p,2]*(U[n,w,2]) );
                    lrn[p,n,w] = (-pow_a[p] + pow_a[p] * pow(ses_count,-pow_b[n,p]));
                }
                noise[1,n,w] = sigma_pr_r[1]*params_smb_pr[n,w,1];
                noise[2,n,w] = sigma_pr_r[2]*params_smb_pr[n,w,2];
                noise[3,n,w] = sigma_pr_r[3]*params_smb_pr[n,w,3];

               
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



/////// generated quantities {
    // For posterior predictive check
///////    real y_pred_all_weeks_smb[Tr_max_smb, W, N];

///////    real log_lik_all_subs[N, W];

///////    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
///////    for (n in 1:N) {
///////        for (w in 1:W) {
///////            for (t in 1:Tr_max_smb) {
///////                y_pred_all_weeks_smb[t, w, n] = -2;
///////            }
///////        }
///////    }

    // here start predictive checks and logliklhd calculations
///////    for (n in 1:N) {
///////        for (w in 1:W) {

///////            log_lik_all_subs[n,w] = 0;

///////            if (idx_smb_obs[n,w] != 0) {                                // if week exists
///////                    matrix[Tr_smb[n,w], P_smb] X_smb;
///////                    X_smb = to_matrix(x_smb[n,w,:Tr_smb[n,w],]);
///////                    y_pred_all_weeks_smb[:Tr_smb[n,w], w, n] = bernoulli_rng(Phi(X_smb * to_vector(params_smb[n,w,])));   
///////                    log_lik_all_subs[n,w] += bernoulli_lpmf(y_smb[n,w,:Tr_smb[n,w]] | Phi(X_smb * to_vector(params_smb[n,w,])));               
///////            }                         
///////        }
///////    }


     // Here we save the learning curves per week
///////    real learning_curve_w_V_smb[N,W]; // Value of the mean
///////    real learning_curve_w_signVoverTU_smb[N,W]; // Value of the mean
///////    real learning_curve_w_RU_smb[N,W]; // Value of the mean

///////    for (n in 1:N) {
///////        real ses_count = 2.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
///////        for (w in 1:W) {
///////            learning_curve_w_V_smb[n,w] = -2;
///////            learning_curve_w_signVoverTU_smb[n,w] = -2;
///////            learning_curve_w_RU_smb[n,w] = -2;
///////           if (idx_smb_obs[n,w] != 0) {
///////                learning_curve_w_V_smb[n,w] =           (mu_pr_sub[1,n] + pow_a[1] * pow(ses_count,-pow_b[n,1]));     // smb
///////                learning_curve_w_signVoverTU_smb[n,w] = (mu_pr_sub[2,n] + pow_a[2] * pow(ses_count,-pow_b[n,2]));    // smb
///////                learning_curve_w_RU_smb[n,w] =          (mu_pr_sub[3,n] + pow_a[3] * pow(ses_count,-pow_b[n,3]));     // smb
///////                ses_count = ses_count + 1.0;
///////            }    
///////        }
///////    }

/////// }