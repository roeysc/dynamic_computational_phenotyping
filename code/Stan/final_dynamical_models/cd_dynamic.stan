// In this model, the priors are N(0,1), that is, not based on the literature.



functions {

    real partial_sum(int[,,] choice_itc, int start, int end,

			        int[,,] choice_cd, real[,] criterion_cd, real[,] sigma_cd, real[,] criterion_slope_cd, real[,] sigma_slope_cd, 
			        int[,] idx_cd_obs, int[,,] Nb_cd, int[,,] Tar_cd,  real[,,] delta_cd, int[,]  Tr_cd, int W) {

	real lt = 0;
	for (n in 1:(end-start+1)) {

		int s = start + (n - 1); 

		for (w in 1:W) {

 	        if (idx_cd_obs[s,w] != 0) { 
                for (t in 1:Tr_cd[s,w]) {
                    // Figure out what block this is based on the number of targets
        	        real factor = Nb_cd[s,w,t]-3; // // Factor is the factor by which we need to multiply sigma and critrion, which are assumed to increase linearly with the number of items. The baseline is 3, which is the minimal number of items in the task.

                    real f = 0.5 * (erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor))/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))) + erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor))/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))));       

                    // generate posterior prediction for current trial
                    if (delta_cd[s,w,t]==0){
                	    real h = 0.5 * (erfc((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))) + erfc((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))));
                        real p = pow((1-h), Tar_cd[s,w,t]) * pow((1-f),(Nb_cd[s,w,t]-Tar_cd[s,w,t]));   
                        lt += bernoulli_lpmf(choice_cd[s,w,t] | 1.0-p); // Becuase choice_cd is what's being sliced over, it's indexed with n and not s
                	} else {
                	    real h = 0.5 * (erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)-1)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))) + erfc(((criterion_cd[s,w]+criterion_slope_cd[s,w]*factor)+1)/(sqrt(2)*(sigma_cd[s,w]+sigma_slope_cd[s,w]*factor))));
                        real p = pow((1-h), Tar_cd[s,w,t]) * pow((1-f),(Nb_cd[s,w,t]-Tar_cd[s,w,t]));   
                        lt += bernoulli_lpmf(choice_cd[s,w,t] | 1.0-p);
                    }
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
    int pressed_gng[N,W,Bl, Tr_max_gng];
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
    int choice_lt[N, W, Tr_max_lt]; // Risky=0; Safe=1. Roey: Why is -1 an option here?   

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
    int num_par = P_cd;
    int exo_q_num_use = 2; // Number of exo effects to use
}

parameters {
    vector[num_par] mu_pr_raw;
    vector<lower=0>[num_par] sigma_pr;
    real mu_pr_sub[num_par, N];
    real<lower=0> sigma_pr_r[num_par];
    vector[2] pow_a_pr; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training
    vector[2] mu_pow_b_pr; // A prior for pow_b, for each one of the phenotype patameters
    vector<lower=0>[2] sigma_pow_b_pr; // Sigma of pow_b, to allow for between-subject variability for each one of the phenotype patameters
    
    //vector[num_par] w_prv_pr; // Weight for the previous week in the slow-dynamics model (one parameter shared across all subjects)

    real criterion_cd_pr[N,W]; //    4 criteria, one for each block
    real sigma_cd_pr[N,W];     // 4 sigmas, one for each block
    real criterion_slope_cd_pr[N,W]; //    4 criteria, one for each block
    real sigma_slope_cd_pr[N,W];     // 4 sigmas, one for each block
    real pow_b_pr[N,2]; // For modelling param = param_baseline + a*(time)^b (with 'time' being a count of sessions)
    
    matrix[2, exo_q_num_use] A_pr; //EXOMODEL
}

transformed parameters {
    vector[num_par] mu_pr;
    mu_pr[1] = mu_pr_raw[1];
    mu_pr[2] = mu_pr_raw[2];
    mu_pr[3] = mu_pr_raw[3]; // normal(0,1)
    mu_pr[4] = mu_pr_raw[4]; // normal(0,1)
    

    real criterion_cd[N,W]; //    4 criteria, one for each block
    real sigma_cd[N,W];     // 4 sigmas, one for each block
    real criterion_slope_cd[N,W]; //    4 criteria, one for each block
    real sigma_slope_cd[N,W];     // 4 sigmas, one for each block
    real pow_a_all_subjects[2,N];
    matrix[2, exo_q_num_use] A_all_subjects[N]; //EXOMODEL
    real<lower=0> pow_b[N,2]; // For modelling param = param_baseline  + a*(time)^b (with 'time' being a count of sessions). This can be both positive and negative, depending on whether this parameter tends to increase/decrease with training

    //vector<lower=0>[2] w_prv; // Weight for the previous week in the slow-dynamics model (one parameter shared across all subjects)
    
    for (n in 1:N) {
        for (p in 1:2){
            pow_b[n,p] = exp(-2) + (exp(1)-exp(-2))*Phi(mu_pow_b_pr[p] + sigma_pow_b_pr[p]*pow_b_pr[n,p]); // Constrain b to values that don't lead to flat learning curves (which effectively remove the learning from the learning term, causeing identifiability issues)
        }
        real ses_count = 1.0;  // Session count, to be used as "time" in the power law. This version starts from 2.0 so that we already  have an effect for very large exponent
        int first_session; // Boolean
        first_session = 1;
        
         // Create the subject-specific learning term pow_a
        vector[2] pow_a;
        pow_a[1] = (2.0*Phi_approx(pow_a_pr[1])-1.0).*mu_pr_sub[1,n]; // each pow_a_frac is a fraction of the corresponding mu_pr_sub (multiplied by 2 to allow for greater learning)
        pow_a[2] = (2.0*Phi_approx(pow_a_pr[2])-1.0).*mu_pr_sub[3,n]; // each pow_a_frac is a fraction of the corresponding mu_pr_sub (multiplied by 2 to allow for greater learning)
        pow_a_all_subjects[1,n] = pow_a[1];
        pow_a_all_subjects[2,n] = pow_a[2];
        
        // Create the subject-specific exogenous term
        matrix[2, exo_q_num_use] A; //EXOMODEL
        for (q in 1:exo_q_num_use){
            A[1,q] = (2.0*Phi_approx(A_pr[1,q])-1.0).*mu_pr_sub[1,n]; // each exogenous loading is a fraction of the corresponding mu_pr
            A[2,q] = (2.0*Phi_approx(A_pr[2,q])-1.0).*mu_pr_sub[3,n]; // each exogenous loading is a fraction of the corresponding mu_pr
        }
        A_all_subjects[n] = A;
        

        
        real prev_criterion_cd; // Previous value in unconstrained space
        real prev_sigma_cd; 
        
        for (w in 1:W) {
            criterion_cd[n,w] = 0; // initializae, so that we don't have undefined parameters
            criterion_slope_cd[n,w] = 0;
            sigma_cd[n,w] = 0;
            sigma_slope_cd[n,w] = 0;
            if (idx_cd_obs[n,w] != 0) { // if week exists
                //if (first_session == 1) {
                //    prev_criterion_cd       = mu_pr_sub[1,n];
                //    prev_sigma_cd           = mu_pr_sub[3,n];
                //    first_session = 0;
                //}
                criterion_cd[n,w] =       exp(mu_pr_sub[1,n]   +   (A[1,]*to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1]))   +   sigma_pr_r[1]*criterion_cd_pr[n,w]);
                criterion_slope_cd[n,w] = exp(mu_pr_sub[2,n]   +   sigma_pr_r[2]*criterion_slope_cd_pr[n,w]);
                sigma_cd[n,w] =           exp(mu_pr_sub[3,n]   +   (A[2,]*to_vector(U[n,w,1:exo_q_num_use]))   +    (- pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2]))   +   sigma_pr_r[3]*sigma_cd_pr[n,w]);
                sigma_slope_cd[n,w] =     exp(mu_pr_sub[4,n]   +   sigma_pr_r[4]*sigma_slope_cd_pr[n,w]);
                
                //prev_criterion_cd       =     w_mu[1]*(mu_pr_sub[1,n])   +   w_prv[1]*(prev_criterion_cd)         +   w_exo[1]*(mu_pr_sub[1,n] + A[1,]*to_vector(U[n,w,]))   +    w_lrn[1]*(mu_pr_sub[1,n] - pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1]))   +   sigma_pr_r[1]*criterion_cd_pr[n,w];
                //prev_sigma_cd           =     w_mu[2]*(mu_pr_sub[3,n])   +   w_prv[2]*(prev_sigma_cd)             +   w_exo[2]*(mu_pr_sub[3,n] + A[2,]*to_vector(U[n,w,]))   +    w_lrn[2]*(mu_pr_sub[3,n] - pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2]))   +   sigma_pr_r[3]*sigma_cd_pr[n,w];
                ses_count = ses_count+1.0;
            }
        }
    }
}

model {
    // hyper parameters
    mu_pr_raw ~ normal(0, 1.0); // normal(-1.5, 0.51); // So sigma_cd ~ normal(0.26, 0.15). Based on Wilken and Ma (2004), Figure 5. I extracted the values and calculated a linear fit, using the 95% CI to get an estimate of the STD (by taking CI/4).  Then, I used fminsearch to find the values of mu and sigma in the unconstrained space (as mentioned here: https://docs.google.com/document/d/1lB45Ixnj40xqTFhZz17mXz3vAaX8Pv3y1CVZx6HqlXg/edit?usp=sharing)

    sigma_pr ~ normal(0, sqrt(1.0)); // We could shrink this, because we assume subject are overall similar to each other
    pow_a_pr ~ normal(0, 1.0); // We could shrink this, because we want dynamical parameters to be close to zero
    mu_pow_b_pr  ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. We do shrink pow_a (C in the overlead notation), so I think that's enough
    
    //w_prv_pr ~ normal(0, 1.0);
    to_vector(A_pr) ~ normal(0, sqrt(1.0)); // ***  We shrink this, because we want dynamical parameters to be close to zero
    
    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr_raw[p], sigma_pr[p]); // ROEY: USED MU_PR_RAW HERE JUST TO MAKE IT SAMPLE, PROBABLY BEST TO MOVE THIS TO THE TRANSFORMED PARAMETERS OR SOMETHING
    }

    to_vector(sigma_pr_r) ~ normal(0, 1.0); // Here we DON'T shrink the week-to-week noise
    sigma_pow_b_pr ~ normal(0, 1.0); // We DON'T shrink this, becuase we want to allow freedom in the learning dynamics. &&&TODO: MAYBE WE DO WANT TO SHRINK THIS UNDER THE ASSUMPTION THAT SUBJECTS ARE SIMILAR TO EACH OTHER?

    // individual parameters w/ Matt trick
    to_vector(to_matrix(criterion_cd_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(criterion_slope_cd_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(sigma_cd_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(sigma_slope_cd_pr)) ~ normal(0, 1.0);
    // same for learning curves
    to_vector(to_matrix(pow_b_pr)) ~ normal(0, 1.0);

    	target += reduce_sum(partial_sum, choice_itc, 1, choice_cd, criterion_cd, sigma_cd, criterion_slope_cd, 
    	                    sigma_slope_cd, idx_cd_obs, Nb_cd, Tar_cd, delta_cd, Tr_cd, W);
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
        vector[2] pow_a;
        pow_a = to_vector(pow_a_all_subjects[,n]);

        // Get the subject-specific exogenous term
        matrix[2, exo_q_num_use] A; //EXOMODEL
        A = A_all_subjects[n];
        
        for (w in 1:W) {
            if (idx_cd_obs[n,w] != 0) { // if week exists
                // Fill in the terms we will use to calculate the variance explained of each term                
                exo1[1,n,w] = (A[1,1]*(U[n,w,1]) );
                exo2[1,n,w] = (A[1,2]*(U[n,w,2]) );
                lrn[1,n,w] = (-pow_a[1] + pow_a[1] * pow(ses_count,-pow_b[n,1]));
                noise[1,n,w] = sigma_pr_r[1]*criterion_cd_pr[n,w];

                exo1[2,n,w] = (A[2,1]*(U[n,w,1]) );
                exo2[2,n,w] = (A[2,2]*(U[n,w,2]) );
                lrn[2,n,w] = (-pow_a[2] + pow_a[2] * pow(ses_count,-pow_b[n,2]));
                noise[2,n,w] = sigma_pr_r[3]*sigma_cd_pr[n,w]; // Indeed, it's sigma_pr_r[3] (see above in the model)

                ses_count = ses_count+1.0;
            }
        
        }
    }
    
    // *** Get the general exo and learning terms, to calculate their probability of direction
    vector[2] pow_a_frac;
    matrix[2,exo_q_num_use] A_frac;
    real sgn; // we multiply the terms by the sign of the parameter's baseline so we know of the effect is increasing or decreasing the parameter (for leaerning, s negative sign means increasing - checl out the opwerpower law formulation to see why) 

    for (p in 1:2){
        if (mu_pr_raw[p]<0){
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
 