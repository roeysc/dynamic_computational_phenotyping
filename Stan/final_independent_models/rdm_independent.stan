// wiener_rng function was generously provided by Sean Pinkney; see https://discourse.mc-stan.org/t/alternative-to-wiener-rng-in-stan/28443


functions {
// The following functions are required for the RDM task (for sampling from the Wiener distribution)
real fs_cdf(real t, real a) {
    if (a < 1) {
      reject("a must be >= 1, found a = ", a);
    }

    return erfc(inv_sqrt(2 * a * t));
  }

  vector make_vars(real mu) {
    real mu2 = pow(mu, 2);
    real t_tilde = 0.12 + 0.5 * exp(-mu2 / 3);
    real a = (3 + sqrt(9 + 4 * mu2)) / 6;
    real sqrtamu = sqrt((a - 1) * mu2 / a);
    real fourmu2pi = (4 * mu2 + pi() ^ 2) / 8;
    real Cf1s = sqrt(a) * exp(-sqrtamu);
    real Cf1l = pi() / (4 * fourmu2pi);
    real CF1st = Cf1s * fs_cdf(t_tilde | a);
    real F1lt = -expm1(-t_tilde * fourmu2pi);
    real F1inf = CF1st + Cf1l * (1 - F1lt);

    return [mu2, //.......1
            t_tilde, //.. 2
            a, //.........3
            sqrtamu, //...4
            fourmu2pi, //.5
            Cf1s, //......6
            Cf1l, //......7
            CF1st, //.....8
            F1lt, //......9
            F1inf]'; //...10
  }

  int acceptt_rng(real t_star, real ft, real c) {
    if (c <= 0.06385320297074884) {
      reject("c is ", c);
    }
    if (is_nan(c)) {
      reject("c is nan!");
    }
    real z = ft * uniform_rng(0, 1);
    real b = exp(-c);
    int k = 3;

    while (1) {
      if (z > b) {
        return 0;
      }
      b -= k * exp(-c * k ^ 2);
      if (z < b) {
        return 1;
      }
      k += 2;
      b += k * exp(-c * k ^ 2);
      k += 2;
    }

    return 0;
  }

  real sample_small_mu_rng(vector vars) {
    real t_star;
    real pi_sq = pi() ^ 2;

    real mu2 = vars[1];
    real a = vars[3];
    real sqrtamu = vars[4];
    real fourmu2pi = vars[5];
    real Cf1s = vars[6];
    real Cf1l = vars[7];
    real CF1st = vars[8];
    real F1lt = vars[9];
    real F1inf = vars[10];

    int counter_outer = 0;
    while (1) {
      real p = F1inf * uniform_rng(0, 1);

      if (p <= CF1st) {
        t_star = 1. / (2 * a * pow(inv_erfc(p / Cf1s), 2));
        while (0.5 * t_star <= 0.06385320297074884) {
          p = uniform_rng(0.06385320297074884, CF1st);
          t_star = 1. / (2 * a * pow(inv_erfc(p / Cf1s), 2));
        }
        real ft = exp(-1. / (2 * a * t_star) - sqrtamu + mu2 * t_star);
        if (acceptt_rng(t_star, ft, 0.5 * t_star) == 1) {
          return t_star;
        }
      } else {
        t_star = -log1p(-(p - CF1st) / Cf1l - F1lt) / fourmu2pi;
        real pisqt = pi_sq * t_star / 8;
        while (pisqt <= 0.06385320297074884) {
          p = uniform_rng(CF1st, F1inf);
          t_star = -log1p(-(p - CF1st) / Cf1l - F1lt) / fourmu2pi;
          pisqt = pi_sq * t_star / 8;
        }
        if (acceptt_rng(t_star, exp(-pisqt), pisqt) == 1) {
          return t_star;
        }
      }
    }
    return 0;
  }

  real inverse_gaussian_rng(real mu, real mu_sq) {
    real v = pow(std_normal_rng(), 2);
    real z = uniform_rng(0, 1);
    real x = mu + 0.5 * mu_sq * v - 0.5 * mu * sqrt(4 * mu * v + mu_sq * v ^ 2);
    if (z <= (mu / (mu + x))) {
      return x;
    } else {
      return mu_sq / x;
    }
  }

  real sample_large_mu_rng(vector vars) {
    real mu2 = vars[1];
    real t_tilde = vars[2];
    real a = vars[3];
    real sqrtamu = vars[4];
    real fourmu2pi = vars[5];
    real Cf1s = vars[6];
    real Cf1l = vars[7];
    real CF1st = vars[8];
    real F1lt = vars[9];
    real F1inf = vars[10];

    real invabsmu = inv_sqrt(mu2);

    if (t_tilde >= 0.63662) {
      Cf1l = -log(pi() * 0.25) - 0.5 * log(2 * pi());
      Cf1s = 0;
    } else {
      Cf1l = -pi() ^ 2 * t_tilde / 8 + (3. / 2.) * log(t_tilde) + 0.5 * inv(t_tilde);
      Cf1s = Cf1l + 0.5 * log(2 * pi()) + log(pi() * 0.25);
    }

    while (1) {
      real t_star = inverse_gaussian_rng(invabsmu, inv(mu2));
      if (is_nan(t_star)) {
        reject("t_star is nan! ", mu2);
      }
      real one2t = 0.5 * inv(t_star);
      if (t_star <= 2.5) {
        real expone2t = exp(Cf1s - one2t);
        if (expone2t == 0) {
          expone2t = 1e-8;
        }
        if (acceptt_rng(t_star, expone2t, one2t) == 0 || invabsmu < 0.000666) {
          return t_star;
        }
      } else {
        real expone2t = exp(-log(pi() / 4) - 0.5 * log(2 * pi()) - one2t - (3. / 2.) * log(t_star));
        if (acceptt_rng(t_star, expone2t, pi() ^ 2 * t_star / 8) == 0) {
          return t_star;
        }
      }
    }
    return 0;
  }

  real fast_pt_rng(real alpha, real tau, real beta, real delta) {
    real absmu = fabs(delta) ;
    vector[10] vars = make_vars(absmu);
    real pt;

    if (absmu < 1) {
      pt = sample_small_mu_rng(vars);
    } else {
      pt = sample_large_mu_rng(vars);
    }

    return pt;
  }

  vector wiener_rng(real alpha, real tau, real beta, real delta) {
    real t = 0 ;
    real sign_delta = delta > 0 ? 1 : -1;
    real x =  beta * alpha ;
    real mu = fabs(delta);
    real hit_bound;
    vector[2] out;
    int counter = 0;

    if (beta == 0 || beta == 1) {
      return [tau, beta]';
    }

    while (1) {
      real mutheta;
      real xlo =  x ;
      real xhi = alpha - x ;
      // lower bound is 0 
      // upper bound is alpha in stan parmeterization

      // symmetric case, [x - xup, x + xup]
      if (fabs(xlo - xhi) < 1e-6) {
        mutheta = xhi * mu;
        real pt = fast_pt_rng(alpha, tau, beta, xhi * fabs(delta));
        hit_bound = sign_delta == 1 ? inv_logit( 2 * mutheta ) : 1 - inv_logit( 2 * mutheta );
        real bound = uniform_rng(0, 1) < hit_bound ? 1 : 0;
        return [ tau + t  + ( square(xhi)  * pt), bound]';
      // x is closer to upper bound, [x - xup, x + xup]
      } else if (xlo > xhi) {
        mutheta = xhi * mu;
        t += ( square(xhi ) * fast_pt_rng(alpha, tau, beta, xhi* fabs(delta)))  ;
        hit_bound = sign_delta == 1 ? inv_logit( 2 * mutheta ) : 1 - inv_logit( 2 * mutheta );
        if (uniform_rng(0, 1) < hit_bound ) {
          return [tau + t, 1]';
        }
        x -= xhi ;
      } else {
       // x is closer to lower bound, [x - xlo, x + xlo]
        mutheta = xlo * mu ;
        t +=  ( square(xlo ) * fast_pt_rng(alpha, tau,  beta, xlo* fabs(delta) )) ;
        hit_bound = sign_delta == 1 ? inv_logit( 2 * mutheta ) : 1 - inv_logit( 2 * mutheta );
        if (uniform_rng(0, 1) > hit_bound) {
          out[1] = tau + t;
          out[2] = 0 ;
          break;
        }
        x += xlo ;
      }
    }
    return out;
  }


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
    int<upper=1> choice_itc[N,W,Tr_max_itc];     // choice itc - 0 for instant reward, 1 for delayed reward - subj x weeks x trials

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
}

parameters {
    vector[num_par] mu_pr;
    real mu_pr_sub[num_par, N];
    vector<lower=0>[num_par] sigma_pr;
    real<lower=0> sigma_pr_r[num_par];

    real alpha_rdm_pr[N,W];              // rdm
    real delta_rdm_pr[N,W];              // rdm
    real tau_rdm_pr[N,W];                // rdm

}

transformed parameters {
    real tau_rdm[N,W];
    real<lower=0> alpha_rdm[N,W];
    real delta_rdm[N,W];

    
    for (n in 1:N) {
        for (w in 1:W) {
            tau_rdm[n,w] = inv_logit(mu_pr_sub[1,n] + sigma_pr_r[1]*tau_rdm_pr[n,w]) * (minRT_rdm[n,w] - RTbound_rdm - 0.0001) + RTbound_rdm; ;                // rdm
            alpha_rdm[n,w] = exp(mu_pr_sub[2,n] + sigma_pr_r[2]*alpha_rdm_pr[n,w]);          // rdm
            delta_rdm[n,w] = (mu_pr_sub[3,n] + sigma_pr_r[3]*delta_rdm_pr[n,w]);          // rdm
        }
    }

}

model {
    // hyper parameters
    mu_pr[1]  ~ normal(0, 1.0); // We use the value from literature
    mu_pr[2]  ~ normal(0, 1.0); // We use the value from literature
    mu_pr[3]  ~ normal(0, 1.0); // We use the value from literature, but increased the STD by x100, because the choice of minRT to use was quite arbitrary (we used the 5th percentile across the population, but in practice we use the minRT of each session separtely)
    sigma_pr ~ normal(0, 1);

    // subject level parameters
    for (p in 1:num_par) {
            to_vector(mu_pr_sub[p,]) ~ normal(mu_pr[p], sigma_pr[p]);
    }

    to_vector(sigma_pr_r) ~ normal(0, 1);
    
    to_vector(to_matrix(tau_rdm_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(alpha_rdm_pr)) ~ normal(0, 1.0);
    to_vector(to_matrix(delta_rdm_pr)) ~ normal(0, 1.0);

    	target += reduce_sum(partial_sum, choice_itc, 1, RTu_rdm, RTl_rdm, Cohu_rdm, Cohl_rdm, delta_rdm,alpha_rdm, 
    	                    tau_rdm, idx_rdm_obs, Nu_rdm, Nl_rdm, W);
}    

 generated quantities {
    // For posterior predictive check
    real y_pred_all_weeks_rdm[Nu_max_rdm + Nl_max_rdm, W,N];
    real trial_type[Nu_max_rdm + Nl_max_rdm,W, N]; // For RDM, to keep track of correct and incorrect responses in the simulated data

    real log_lik_all_subs[N, W];

    int b;

    // Set all posterior predictions to -2 (avoids NULL values)
    for (n in 1:N) {
        for (w in 1:W) {
            for (t in 1:Nu_max_rdm+Nl_max_rdm) {
                y_pred_all_weeks_rdm[t,w,n] = -2;
                trial_type[t,w,n] = -2;
            }
        }
    }

    // here start predictive checks and logliklhd calculations
    for (n in 1:N) {
        for (w in 1:W) {

            log_lik_all_subs[n,w] = 0;            

            if (idx_rdm_obs[n,w] != 0) {

                vector[Nu_rdm[n,w]] delta_cohu = delta_rdm[n,w]*to_vector(Cohu_rdm[n,w,:Nu_rdm[n,w]]);
                vector[Nl_rdm[n,w]] delta_cohl = delta_rdm[n,w]*to_vector(Cohl_rdm[n,w,:Nl_rdm[n,w]]);

                log_lik_all_subs[n,w] += wiener_lpdf(RTu_rdm[n,w,:Nu_rdm[n,w]] | alpha_rdm[n,w], tau_rdm[n,w], 0.5, delta_cohu);
                log_lik_all_subs[n,w] += wiener_lpdf(RTl_rdm[n,w,:Nl_rdm[n,w]] | alpha_rdm[n,w], tau_rdm[n,w], 0.5, -delta_cohl); 

                //generate posterior predictions for current trials
                for (t in 1:Nu_rdm[n,w]) {
        	        vector[2] out = wiener_rng(alpha_rdm[n,w], tau_rdm[n,w], 0.5, delta_cohu[t]);
                    y_pred_all_weeks_rdm[t,w,n] = out[1]; // predicted RT
                    trial_type[t,w,n] = out[2]; // correct or incorrect response
	            }

                for (t in 1:Nl_rdm[n,w]) {
            	    vector[2] out = wiener_rng(alpha_rdm[n,w], tau_rdm[n,w], 0.5,delta_cohl[t]);
                    y_pred_all_weeks_rdm[Nu_rdm[n,w] + t,w,n] = out[1]; // predicted RT
                    trial_type[Nu_rdm[n,w] + t,w,n] = out[2];  // correct or incorrect response
                }
            }    

        }
    }
 }