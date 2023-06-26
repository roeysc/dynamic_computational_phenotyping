function [pheno,parameters_list,pheno_non_zscore,pheno_with_outliers] = comppheno_get_pheno(remove_outliers,outliers_output_text_file,model,gng_model)
% comppheno_get_pheno_hierarchical returns a matrix of phenotype parameters for all sessions,
% after applying transformations (e.g., log() for k_itc) and z-score.
% The output is arranged as parameter x subjects x week
%
% Missing values are replaced with NaN
% Outlier removal is not recommended. However, if outliers are removed,
% they are replaced with NaN.
% Outlier removal might be necessary when using parameter estimation
% methods that are not hierarchical, like session-wise MLE. In these cases,
% there in no natural pooling of parameter values towards the
% participant-level mean (and the population-level mean), and outliers are
% more likely to occur.
%
% 'remove_outliers' controls the definition of outliers to be removed:
% 1 = Outliers are elements more than 3 standardizes MAD from the median
% 2 = Outliers are elements more than 1.5 IQR from the median
% 3 = Outliers are the first and last percentile
% 4 = No outlier removal
% 'outliers_output_text_file' is optional, for saving how many outliers were
% found in the data
%
% 'model' controls the statistical model used to fit the computional
% parameters:
%
% 1. 'independent' - the baseline independent model, where each
% session is sampled from a participant-specific distribution, whose
% parameters are fit using higher population-level parameters.
%
% 2. 'dynamic' - the dynamic statistical model, which has assumptions about the
% structured variation of the phenotype across time (practice effects,
% affective state effects etc.)
%
% 3. 'simulated_data_based_on_phenotype_of_first_session' - this is used for
% calculating the upper bound of the ICC values. Under the ground-truth
% conditions of a stable phenotype within participant (used to simulate the
% data of all sessions), it allows to test whether ICC is high (meaning
% that there are enough trials in the experimental paradigm and that the
% fitting procedure is robust enough), such that the recovered phenotype is
% indeed stable like the ground-truth one.
%
% 'gng_model' controls which model for the Go/NoGo task will be used:
% 1. 'modified' = The revised model we propose. In this model, there are
%      two 'rho' parameters. Ine is the effective size of rewards/punishments,
%      and one is the effective size of neutral outcome.
% 2. 'original' = The original Guitart-Masip (2012) model.

%% Set directories
comppheno_set_dirs % Load the comppheno_analysis_dir variable

%% Set default arguments
if ~exist('model','var')
    model = 'independent';
end

if ~exist('gng_model','var')
    gng_model = 'modified';
end

exclude_gng_subejcts_with_accuracy_below_55 = true; % These are the non-learners (throughout all 12 sessions), akin to the non-learners defined in Guitart-Masip et al. (2012)
exclude_lt_subjects_who_choose_safe_more_than_80_percent = true; % Remove subjects

%% Set input directory
if strcmp(model,'independent')
    phenotype_dir = fullfile(comppheno_analysis_dir,'phenotype_independent_model');
elseif strcmp(model,'dynamic')
    phenotype_dir = fullfile(comppheno_analysis_dir,'phenotype_dynamic_model');
elseif strcmp(model,'simulated_data_based_on_phenotype_of_first_session')
    phenotype_dir = fullfile(comppheno_analysis_dir,'phenotype_simulated_data_based_on_phenotype_of_first_session');
elseif strcmp(model,'simulated_data_based_on_the_full_phenotype')
    phenotype_dir = fullfile(comppheno_analysis_dir,'phenotype_simulated_data_based_on_the_full_phenotype');
end

%% Set parameter albels
params_cd = {'criterion_cd','criterion_slope_cd','sigma_cd','sigma_slope_cd'};
params_itc = {'itc_k','itc_beta'};
params_lt = {'risk_lt','beta_lt'};
params_nc = {'weber_nc'};
params_smb = {'params_smb_1','params_smb_2','params_smb_3'}; % Matlab changes variables names from having . to having _, so we use _ here
params_rdm = {'alpha_rdm','delta_rdm','tau_rdm'};
if strcmp(gng_model,'modified')
    params_gng = {'b_gng','pi_gng','xi_gng','ep_gng','rho_rew_pun_gng','rho_neut_gng'};
elseif strcmp(gng_model,'original')
    params_gng = {'b_gng','pi_gng','xi_gng','ep_gng','rho_gng'};
end
subjects = comppheno_get_subjects();

%% Files and directory settings
% These are the csv phenotype files of the fitted Stan models. They are
% created by openFile_roey.py after running the models in Stan.
nc_file = fullfile(phenotype_dir,sprintf('nc_parameters_%s.csv',model));
cd_file = fullfile(phenotype_dir,sprintf('cd_parameters_%s.csv',model));
itc_file = fullfile(phenotype_dir,sprintf('itc_parameters_%s.csv',model));
lt_file = fullfile(phenotype_dir,sprintf('lt_parameters_%s.csv',model));
smb_file = fullfile(phenotype_dir,sprintf('smb_parameters_%s.csv',model));
gng_file = fullfile(phenotype_dir,sprintf('gng_parameters_%s.csv',model));
rdm_file = fullfile(phenotype_dir,sprintf('rdm_parameters_%s.csv',model));


%% Fill in all parameters task by task
pI = 1;
parameters_list = {}; % Fill a list of parameters, to keep track of which column in which variable

%% cd gng itc lt nc smb rdm
pheno = [];
parameters_list = {};

missing_sessions = comppheno_get_missing_sessions('gng');

try
    pheno_struct = comppheno_read_hierarchical_csv(gng_file);
catch
    warning('Failed to read the phenotype for GNG.')
    % If file doesn't exist, fill with NaNs
    for pI = 1:length(params_gng)
        pheno_struct.(params_gng{pI}) = nan(size(missing_sessions));
    end
end


if exclude_gng_subejcts_with_accuracy_below_55 % Arrange the 66 remaining subjects in a 90-long structure
    idx = 1:90;
    idx_excluded = [1, 2, 3, 6, 8, 9, 15, 16, 17, 28, 33, 34, 37, 39, 44, 59, 60, 68, 70, 71, 74, 76, 78, 86]; % These are the subjects we excluded during the model fit
    idx(ismember(idx,idx_excluded))=[];

    for pI = 1:length(params_gng)
        pheno_struct_new.(params_gng{pI}) = nan(size(missing_sessions));
        for sI = 1:length(idx)
            pheno_struct_new.(params_gng{pI})(idx(sI),:) = pheno_struct.(params_gng{pI})(sI,:);
        end
    end
    pheno_struct = pheno_struct_new;
end


for pI = 1:length(params_gng)
    tmp = pheno_struct.(params_gng{pI});
    tmp(missing_sessions) = nan;
    pheno(end+1,:,:) = tmp;
    parameters_list{end+1} = params_gng{pI};
end
% return % &&&

missing_sessions = comppheno_get_missing_sessions('cd');
try
    pheno_struct = comp_pheno_read_hierarchical_csv(cd_file);
catch
    warning('Failed to read the phenotype for CD.')
    % If file doesn't exist, fill with NaNs
    for pI = 1:length(params_cd)
        pheno_struct.(params_cd{pI}) = nan(size(missing_sessions));
    end
end
for pI = 1:length(params_cd)
    tmp = pheno_struct.(params_cd{pI});
    tmp(missing_sessions) = nan;
    pheno(end+1,:,:) = tmp;
    parameters_list{end+1} = params_cd{pI};
end

missing_sessions = comppheno_get_missing_sessions('itc');
try
    pheno_struct = comp_pheno_read_hierarchical_csv(itc_file);
catch
    warning('Failed to read the phenotype for ITC.')
    % If file doesn't exist, fill with NaNs
    for pI = 1:length(params_itc)
        pheno_struct.(params_itc{pI}) = nan(size(missing_sessions));
    end
end
for pI = 1:length(params_itc)
    tmp = pheno_struct.(params_itc{pI});
    tmp(missing_sessions) = nan;
    pheno(end+1,:,:) = tmp;
    parameters_list{end+1} = params_itc{pI};
end

missing_sessions = comppheno_get_missing_sessions('lt');
try
    pheno_struct = comp_pheno_read_hierarchical_csv(lt_file);
catch
    warning('Failed to read the phenotype for LT.')
    % If file doesn't exist, fill with NaNs
    for pI = 1:length(params_lt)
        pheno_struct.(params_lt{pI}) = nan(size(missing_sessions));
    end
end

if exclude_lt_subjects_who_choose_safe_more_than_80_percent  % Arrange the 66 remaining subjects in a 90-long structur
    idx = 1:90;
    idx_excluded = [1,4,7,8,9,11,14,15,19,23,24,25,28,32,35,36,37,47,50,53,58,60,61,65,66,68,70,74,77,83,84,88];
    idx(ismember(idx,idx_excluded))=[];
    for pI = 1:length(params_lt)
        pheno_struct_new.(params_lt{pI}) = nan(size(missing_sessions));
        for sI = 1:length(idx)
            pheno_struct_new.(params_lt{pI})(idx(sI),:) = pheno_struct.(params_lt{pI})(sI,:);
        end
    end
    pheno_struct = pheno_struct_new;
end

for pI = 1:length(params_lt)
    tmp = pheno_struct.(params_lt{pI});
    tmp(missing_sessions) = nan;
    pheno(end+1,:,:) = tmp;
    parameters_list{end+1} = params_lt{pI};
end

missing_sessions = comppheno_get_missing_sessions('nc');
try
    pheno_struct = comp_pheno_read_hierarchical_csv(nc_file);
catch
    % If file doesn't exist, fill with NaNs
    for pI = 1:length(params_nc)
        pheno_struct.(params_nc{pI}) = nan(size(missing_sessions));
    end
end
for pI = 1:length(params_nc)
    tmp = pheno_struct.(params_nc{pI});
    tmp(missing_sessions) = nan;
    pheno(end+1,:,:) = tmp;
    parameters_list{end+1} = params_nc{pI};
end

missing_sessions = comppheno_get_missing_sessions('smb');
try
    pheno_struct = comp_pheno_read_hierarchical_csv(smb_file);
catch
    warning('Failed to read the phenotype for smb.')
    % If file doesn't exist, fill with NaNs
    pheno_struct.params_smb = nan(size(missing_sessions,1),size(missing_sessions,2),3);
end

% smb parameter 1
tmp = squeeze(pheno_struct.params_smb(:,:,1));
tmp(missing_sessions) = nan;
pheno(end+1,:,:) = tmp;
parameters_list{end+1} = 'w_V_smb';
% smb parameter 2
tmp = squeeze(pheno_struct.params_smb(:,:,2));
tmp(missing_sessions) = nan;
pheno(end+1,:,:) = tmp;
parameters_list{end+1} = 'w_signVoverTU_smb';
% smb parameter 3
tmp = squeeze(pheno_struct.params_smb(:,:,3));
tmp(missing_sessions) = nan;
pheno(end+1,:,:) = tmp;
parameters_list{end+1} = 'w_RU_smb';


missing_sessions = comppheno_get_missing_sessions('rdm');
try
    pheno_struct = comp_pheno_read_hierarchical_csv(rdm_file);
catch
    warning('Failed to read the phenotype for RDM.')
    % If file doesn't exist, fill with NaNs
    for pI = 1:length(params_rdm)
        pheno_struct.(params_rdm{pI}) = nan(size(missing_sessions));
    end
end
for pI = 1:length(params_rdm)
    tmp = pheno_struct.(params_rdm{pI});
    tmp(missing_sessions) = nan;
    pheno(end+1,:,:) = tmp;
    parameters_list{end+1} = params_rdm{pI};
end


%% NC (Numerosity comparison)
% Subject 12 has very low accuracy in the NC task across all sessions, and we remove them here
idx = find(~cellfun(@isempty,strfind(parameters_list,'weber_nc')));
pheno(idx,12,:) = nan;

%% Parameter transformations
warning('Transforming ITC''s k parameter to log(k)');
idx = find(strcmp(parameters_list,'itc_k'));
pheno(idx,:,:) = log(pheno(idx,:,:));

%% Remove outliers
fid = fopen(outliers_output_text_file,'w+');
pheno_with_outliers = pheno; % Return this as well for quality checks
if remove_outliers == 1
    disp('Removing phenotype outliers using 3 standardized MAD of the median')
    fprintf(fid,'Removing phenotype outliers using 3 standardized MAD of the median\n');
    for pI = 1:size(pheno,1)
        if ismember(parameters_list{pI},{'itc_beta','beta_lt'})
            disp(['Skipping outlier removal for ' parameters_list{pI} ', which is tightly bounded and also bi-modal']);
            fprintf(fid,'Skipping outlier removal for risk_lt, which is tightly bounded\n');
            continue
        end
        if ismember(parameters_list{pI},{'risk_lt','pi_gng','xi_gng'})
            disp(['Skipping outlier removal for ' parameters_list{pI} ', which is tightly bounded']);
            fprintf(fid,'Skipping outlier removal for risk_lt, which is tightly bounded\n');
            continue
        end
        tmp = pheno(pI,:,:);
        tmp = tmp(:);
        outlier_flag = isoutlier(tmp,'median','ThresholdFactor',3);
        tmp(outlier_flag) = nan;
        tmp = reshape(tmp,size(pheno,2),size(pheno,3));
        pheno(pI,:,:) = tmp;

        disp(sprintf('Found %g outliers (%.2g prcnt) for %s',sum(outlier_flag), 100*sum(outlier_flag)/numel(outlier_flag), parameters_list{pI}))
        fprintf(fid,sprintf('Found %g outliers (%.2g prcnt) for %s\n',sum(outlier_flag), 100*sum(outlier_flag)/numel(outlier_flag), parameters_list{pI}));

    end
end

if remove_outliers == 2
    disp('Removing phenotype outliers using 1.5 IQR for each parameter')
    fprintf(fid,'Removing phenotype outliers using 1.5 IQR for each parameter\n');
    for pI = 1:size(pheno,1)
        if strcmp(parameters_list{pI},'risk_lt')
            disp('Skipping outlier removal for risk_lt, which is tightly bounded');
            fprintf(fid,'Skipping outlier removal for risk_lt, which is tightly bounded\n');
            continue
        end
        tmp = pheno(pI,:,:);
        tmp = tmp(:);
        outlier_flag = isoutlier(tmp,'quartiles');
        tmp(outlier_flag) = nan;
        tmp = reshape(tmp,size(pheno,2),size(pheno,3));
        pheno(pI,:,:) = tmp;

        disp(sprintf('Found %g outliers (%.2g prcnt) for %s',sum(outlier_flag), 100*sum(outlier_flag)/numel(outlier_flag), parameters_list{pI}))
        fprintf(fid,sprintf('Found %g outliers (%.2g prcnt) for %s\n',sum(outlier_flag), 100*sum(outlier_flag)/numel(outlier_flag), parameters_list{pI}));
    end
end

% Methd 3 for removing outliers is below

%% Also return a z-scored version of each parameter separately across all subjects and sessions
pheno_non_zscore = pheno; % Return this as well for quality checks
for pI = 1:size(pheno,1)
    x = pheno(pI,:,:);
    x = x(:);
    x_zscore = (x - nanmean(x))/nanstd(x);
    tmp = reshape(x_zscore,size(pheno,2),size(pheno,3));
    pheno(pI,:,:) = tmp;
end

%% Remove outliers using method 3
% Remove the outliers based on 1st and 99th percentiles (after z-scoring)
if remove_outliers == 3
    prctile_bounds = prctile(pheno(:),[1,99]);
    pheno_interpolated(pheno<prctile_bounds(1))=prctile_bounds(1);
    pheno_interpolated(pheno>prctile_bounds(2))=prctile_bounds(2);
    fprintf(fid,'Removing phenotype outliers by replacing etreme values with the 1st and 99th percentile across all sessions and parameters (after z-scoring)')
end

fclose(fid)