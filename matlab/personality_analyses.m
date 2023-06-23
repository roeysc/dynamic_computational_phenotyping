addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/publish_code/matlab')

%% load personality data

correlate_with_phenotype_mean_or_std = 'baseline'; % mean, median or std or TotalDifference for the difference between week 12 and week 1. 'baseline' (first sessoin)
factor_num = 6; % Number of personality factors to use (3 or 6)
sub_num = 90;

%% Load data
fname = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/pesonality_Barratt_all_subs.csv';
%  fname = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/pesonality_RiskTaking_all_subs.csv';

data_orig = readtable(fname);
subjects = comp_pheno_get_subjects();

% Separate data to in-sample (the final people included in the study) and
% out of sample
idx_in_sample = [];
for sI = 1:length(subjects)
    sub = subjects{sI};
    idx = find(strcmp(data_orig.subjectId,sub));
    data(sI,:) = data_orig(idx,:);
    idx_in_sample(end+1) = idx;
end

data_out_of_sample = data_orig;
data_out_of_sample(idx_in_sample,:) = [];


data_mat = table2array(data(:,4:end));


%% Calculated the sub-scales based on the Factor Structure
% See here: https://www.impulsivity.org/measurement/bis11/
if factor_num==6
    factor_names = {'attention','cognitive instability','motor','perseverance','self control','cognitive complexity'};
    factor_q =   {[5 9 11 20 28],[6 24 26],[2 3 4 17 19 22 25],[16 21 23 30],[1 7 8 12 13 14],  [10 15 18 27 29]}; % Question number
    factor_sgn = {[1 -1 1 -1 1], [1 1 1],  [1 1 1 1 1 1 1],    [1 1 1 -1],   [-1 -1 -1 -1 -1 1],[-1 -1 1 1 -1]};  % Sign of each question
elseif factor_num==3
    factor_names = {'attention','motor','non-planning'};
    factor_q =   {[5 9 11 20 28 6 24 26],[2 3 4 17 19 22 25 16 21 23 30],[1 7 8 12 13 14 10 15 18 27 29]}; % Question number
    factor_sgn = {[1 -1 1 -1 1 1 1 1],  [1 1 1 1 1 1 1 1 1 1 -1],   [-1 -1 -1 -1 -1 1 -1 -1 1 1 -1]};  % Sign of each question
end

data_factor = [];
data_mat = table2array(data(:,4:end));

data_mat = data_mat + 1;

for fI = 1:length(factor_names)
    col_idx = factor_q{fI}; % Indices of columns to use
    tmp = data_mat(:,col_idx); % Take only the questions for this factor
    
    for rev = 1:length(factor_sgn{fI})
        if factor_sgn{fI}(rev) == -1
            tmp(:,rev) = 5 - tmp(:,rev);
        end
    end
    data_factor(:,fI) = sum(tmp,2); % Sum the signed responses
end

%% load relative contribution
in_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/4_hierarchical/R_files_for_plots/models_with_no_weights_and_effects_around_zero/exo2';

fig_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs/relative_STD_d';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
if ~exist(fullfile(fig_dir,'pdf'),'dir')
    mkdir(fullfile(fig_dir,'pdf'));
end
export_figures = 0;
 write_all_pd = 0;
 add_iqr = 0;
export_figures_whole_task = 0;
export_learning_figures = 0;
examine_exo_variability = false;
tasks = {'cd','nc','itc','lt','smb','gng','rdm'};

data_factor_orig = data_factor;

for tI = 1:length(tasks)% tasks = {'cd','nc','itc','lt','smb','gng','rdm'};
    
    % Set task to visualize
    task = tasks{tI};
    
    % NOTE: The order of the parameters here must match the order within the Stan code.
    params.itc = {'k','beta'};
    params.nc = {'weber'};
    params.lt = {'rho','beta'};
    params.smb = {'wV','wsignV_over_TU','wRU'};
    params.cd = {'criterion','sigma'};
    params.gng= {'xi','ep','b','pi','rho_rew_pun','rho_neut'};
    params.rdm = {'alpha','delta'};%{'alpha','delta','tau'};


    for pI = 1:length(params.(task))
        
        data_factor = data_factor_orig;
        
        file_name = [task,'_',params.(task){pI},'.mat'];
        tmp = load(fullfile(in_dir,file_name));

        if strcmp(task,'lt')
           idx = [1,4,7,8,9,11,14,15,19,23,24,25,28,32,35,36,37,47,50,53,58,60,61,65,66,68,70,74,77,83,84,88];
           data_factor(idx,:) = [];
        elseif strcmp(task,'gng')
           idx = [1, 2, 3, 6, 8, 9, 15, 16, 17, 28, 33, 34, 37, 39, 44, 59, 60, 68, 70, 71, 74, 76, 78, 86];
           data_factor(idx,:) = [];
        end
        
        for f = 1: factor_num
            md = median(data_factor(:,f));
            upper_idx = (data_factor(:,f) > md);
            lower_idx = (data_factor(:,f) <= md);
          
            exo1_up =  tmp.sd_frac_exo1(upper_idx);
            exo1_low = tmp.sd_frac_exo1(lower_idx);
            
            [h,p] = ttest2(exo1_up, exo1_low);
            
            if p < 0.05
                disp(['exo1 - ', task, ', parameter ', params.(task){pI}, ', factor ', num2str(f), ...
                    ', p=', num2str(p), ' upper mean=', num2str(mean(exo1_up)), ' lower mean=',num2str(mean(exo1_low))]);
            end
            
            exo2_up =  tmp.sd_frac_exo2(upper_idx);
            exo2_low = tmp.sd_frac_exo2(lower_idx);
            
            [h,p] = ttest2(exo2_up, exo2_low);
            
            if p < 0.05
                disp(['exo2 - ', task, ', parameter ', params.(task){pI}, ', factor ', num2str(f), ...
                    ', p=', num2str(p), ' upper mean=', num2str(mean(exo2_up)), ' lower mean=',num2str(mean(exo2_low))]);
            end
        end
                        
            all_data = sum(data_factor,2);
            
            md = median(all_data);
            upper_idx = (all_data > md);
            lower_idx = (all_data <= md);
          
            exo1_up =  tmp.sd_frac_exo1(upper_idx);
            exo1_low = tmp.sd_frac_exo1(lower_idx);
            
            [h,p] = ttest2(exo1_up, exo1_low);
            
            if p < 0.05
                disp('==========================================================')
                disp(['exo1 - ', task, ', parameter ', params.(task){pI},  ...
                    ', p=', num2str(p), ' upper mean=', num2str(mean(exo1_up)), ' lower mean=',num2str(mean(exo1_low))]);
            end
            
            exo2_up =  tmp.sd_frac_exo2(upper_idx);
            exo2_low = tmp.sd_frac_exo2(lower_idx);
            
            [h,p] = ttest2(exo2_up, exo2_low);
            
            if p < 0.05
                disp('==========================================================')
                disp(['exo2 - ', task, ', parameter ', params.(task){pI} ,  ...
                    ', p=', num2str(p), ' upper mean=', num2str(mean(exo2_up)), ' lower mean=',num2str(mean(exo2_low))]);            
            end         
            
    end
end
    