# compheno_get_dynamic_terms_and_relative_contribution reads the Stan output 
# files (4 chains of posterior samples, i.e., model fits) and export for each
# phenotype parameter *.mat files (for visualization in MATLAB) of the time
# series of each dynamic term (practice effects, affective valence, affective 
# arousl and random noise; all in terms of posterior means), their standard 
# deviations (SD), and their relative SD.
#
# Make sure to install all the required packages.

library("bayesplot")
library("ggplot2")
library("rstan")
library("hexbin")

library("bayestestR")
library("posterior")
library("dplyr") # For the "pull" function below
library("R.matlab")

# Set the wanted task here
task <- 'gng' # 'nc','itc','lt','gng','rdm','smb' ('smb' is the two-armed bandit task)
n_sub=90  # Number of subjects

stan_output_dir <-  # SET PATH TO STAN OUTPUT DIRECTORY (WHERE THE POSTERIOR CHAINS WERE SAVED).
comppheno_dir <- '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/publish_code' # SET THIS TO THE MAIN PROJECT DIRECTORY
R_output_dir <- paste(comppheno_dir,'/analysis/phenotype_dynamic_model/phenotype_terms',sep="")
# Get the Stan output csv files. This pattern here assumes that 4 chains were saved with a long name format starting with year, month, day, hour etc. Hence the "20", at the beginning (2023...)
csvfiles <- dir(path = stan_output_dir,pattern=paste(task,'_dynamic-20[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_[1-4].csv',sep=""), full.names = TRUE)

# Set parameter names and adjust the number of subjects in case some subjects were excluded
if (task=='nc'){
  params <- list("nc_weber")
}
if (task=='itc'){
  params <- list("itc_k","itc_beta")
}
if (task=='lt'){
  n_sub=58
  params <- list("lt_rho","lt_beta")
}
if (task=='gng'){
  n_sub=66
  params <- list("gng_xi","gng_ep","gng_b","gng_pi","gng_rho_rew_pun","gng_rho_neut")
}
if (task=='rdm'){
  params <- list("rdm_alpha","rdm_delta")
}
if (task=='smb'){
  params <- list("smb_wV","smb_wsignV_over_TU","smb_wRU")
}
if (task=='cd'){
  params <- list("cd_criterion","cd_sigma")
}


csvfiles # Display the csv files we found, to make sure we actually found the files


tmp  <- do.call(rbind,lapply(csvfiles,function(i){read.delim(i, comment.char = '#', sep=",")})) # Works for multiple csv files

draws <- as_draws(tmp) # "posterior" library

for (pI in 1:length(params)){ 
    ## Convert df to array
    n_samples <- dim(draws)[1] # Usually 4000
    lrn = array(numeric(),c(n_sub,12,n_samples))
    exo1 = array(numeric(),c(n_sub,12,n_samples))
    exo2 = array(numeric(),c(n_sub,12,n_samples))
    noise = array(numeric(),c(n_sub,12,n_samples))
    
    #draws[is.nan(draws)]<- -2
    for (s in 1:n_sub){
      for (w in 1:12){
        exo1[s,w,] <- pull(draws[ , paste("exo1.",pI,".",s,".",w,sep="")])
        exo2[s,w,] <- pull(draws[ , paste("exo2.",pI,".",s,".",w,sep="")])
        lrn[s,w,] <- pull(draws[ , paste("lrn.",pI,".",s,".",w,sep="")])
        noise[s,w,] <- pull(draws[ , paste("noise.",pI,".",s,".",w,sep="")])
      }
    }
    
    
    ## 3) Calculate the SD  of each term, and plot the fraction of the sum of SD (in the MEAN OF THE POSTERIOR)
    # The rationale here is to think of each term as a contributed with some amplitude, and then compare the amplitudes of the different terms
    exo1[exo1==-2] <- NA # Affective valence
    exo2[exo2==-2] <- NA # Affective arousal
    lrn[lrn==-2] <- NA # Practice
    noise[noise==-2] <- NA # Random noise
    
    # Average across samples to get the posterior mean
    exo1_mean = array(numeric(),c(n_sub))
    exo2_mean = array(numeric(),c(n_sub))
    lrn_mean = array(numeric(),c(n_sub))
    noise_mean = array(numeric(),c(n_sub))
    exo1_mean <- apply(exo1,c(1,2),mean,na.rm=TRUE)
    exo2_mean <- apply(exo2,c(1,2),mean,na.rm=TRUE)
    lrn_mean <- apply(lrn,c(1,2),mean,na.rm=TRUE)
    noise_mean <- apply(noise,c(1,2),mean,na.rm=TRUE)
    
    # Get the SD of each term
    sd_exo1 = array(numeric(),c(n_sub))
    sd_exo2 = array(numeric(),c(n_sub))
    sd_lrn = array(numeric(),c(n_sub))
    sd_noise = array(numeric(),c(n_sub))
    sum_sd = array(numeric(),c(n_sub)) # Sum of SDs (not sd of the sum)
    sd_prev_mean  <- apply(prev_mean,1,sd,na.rm=TRUE) # For each subject, variance across the 12 weeks
    sd_exo1_mean   <- apply(exo1_mean,1,sd,na.rm=TRUE) # For each subject, variance across the 12 weeks
    sd_exo2_mean   <- apply(exo2_mean,1,sd,na.rm=TRUE) # For each subject, variance across the 12 weeks
    sd_lrn_mean   <- apply(lrn_mean,1,sd,na.rm=TRUE) # For each subject, variance across the 12 weeks
    sd_noise_mean <- apply(noise_mean,1,sd,na.rm=TRUE) # For each subject, variance across the 12 weeks
    sum_sd = sd_prev_mean + sd_exo1_mean + sd_exo2_mean + sd_lrn_mean + sd_noise_mean # If summing variances
    sum_sd = sd_exo1_mean + sd_exo2_mean + sd_lrn_mean + sd_noise_mean # If summing variances
    #sum_var = sqrt(var_prev_mean^2 + var_exo_mean^2 + var_lrn_mean^2 + var_noise_mean^2) # If summing STD's
    
    # Get the fraction of the sum of SDs 
    sd_frac_exo1 = array(numeric(),c(n_sub))
    sd_frac_exo2 = array(numeric(),c(n_sub))
    sd_frac_lrn = array(numeric(),c(n_sub))
    sd_frac_noise = array(numeric(),c(n_sub))
    sd_frac_exo1   <- sd_exo1_mean/sum_sd
    sd_frac_exo2   <- sd_exo2_mean/sum_sd
    sd_frac_lrn   <- sd_lrn_mean/sum_sd
    sd_frac_noise <- sd_noise_mean/sum_sd
    
    
    #boxplot(as.vector(sd_frac_prev),as.vector(sd_frac_lrn),as.vector(sd_frac_exo),as.vector(sd_frac_noise),ylim=c(0,1))
    boxplot(as.vector(sd_frac_noise),as.vector(sd_frac_lrn),as.vector(sd_frac_exo1),as.vector(sd_frac_exo2),ylim=c(0,1))
    title(paste("sd/sum(sd), param",pI))

    ## Extract the general exo and learning terms, to calculate probability of directoin
    pow_a_frac <- pull(draws[ , paste("pow_a_frac.",pI,sep="")]) # This is the "a" part in the power-law (a*mu_s)*t^-b
    A_frac_1 <- pull(draws[ , paste("A_frac.",pI,".1",sep="")]) # This is the "a" part in (a*mu_s)*Valence
    A_frac_2 <- pull(draws[ , paste("A_frac.",pI,".2",sep="")]) # This is the "a" part in (a*mu_s)*Arousal

    
    pd_pow_a_frac <- pd(pow_a_frac)
    pd_A_frac_1 <- pd(A_frac_1)
    pd_A_frac_2 <- pd(A_frac_2)

    # Reflect the sign of the parameter in the probability of direction (pd).
    if (isTRUE(sum(pow_a_frac<0)>(length(pow_a_frac)/2))){
      pd_pow_a_frac <- -1*pd_pow_a_frac # Note that a negative values here actually leads to an increase in the phenotype parameter (See eq. 8 in the paper)
    }
    if (isTRUE(sum(A_frac_1<0)>(length(pd_A_frac_1)/2))){
      pd_A_frac_1 <- -1*pd_A_frac_1
    }
    if (isTRUE(sum(A_frac_2<0)>(length(pd_A_frac_2)/2))){
      pd_A_frac_2 <- -1*pd_A_frac_2
    }
    
    # Write to matlab file to visualize the phenotype and its sub-terms
    writeMat(paste(R_output_dir,params[pI],".mat",sep=""),sd_frac_exo1=sd_frac_exo1,sd_frac_exo2=sd_frac_exo2,sd_frac_lrn=sd_frac_lrn,sd_frac_noise=sd_frac_noise,sd_frac_prev=sd_frac_prev,prev_mean=prev_mean,exo1_mean=exo1_mean,exo2_mean=exo2_mean,lrn_mean=lrn_mean,noise_mean=noise_mean,pd_pow_a_frac=pd_pow_a_frac[1],pd_A_frac_1=pd_A_frac_1[1],pd_A_frac_2=pd_A_frac_2[1])
}

