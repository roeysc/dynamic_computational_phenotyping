addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/roey_code/dynamical_system/mle')
addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/roey_code/pareto')
addpath('/n/home10/rschurr/Matlab_code/external_toolboxes/altmany-export_fig-3e1e29a'); % For exporting figures
clear

%% Settings
matplotlib_colormaps; % Load colormaps
use_hierarchical_phenotype = true;
use_pyddm = false; % Replace the DDM from Stan with pyddm

correlate_with_phenotype_mean_or_std = 'baseline'; % mean, median or std or TotalDifference for the difference between week 12 and week 1. 'baseline' (first sessoin)
factor_num = 6; % Number of personality factors to use (3 or 6)
sub_num = 90;
alpha_statistic = 0.05; % Threshold for statistics (see random permutations and correction for multiple comparisons below)

export_figures = true;

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

%% Convert to matrix

data_mat = table2array(data(:,4:end));
% data = zscore(data);


%% Run PCA
[coeff, score, latent, tsquared, explained, mu] = pca(data_mat);



figure
plot(cumsum(explained),'-O')


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

for fI = 1:length(factor_names)
    col_idx = factor_q{fI}; % Indices of columns to use

    tmp = data_mat(:,col_idx); % Take only the questions for this factor
    tmp = bsxfun(@times,tmp,factor_sgn{fI}); % Multiple the responses by the correct sign
    data_factor(:,fI) = sum(tmp,2); % Sum the signed responses
end


%% Get the phenotype
if use_hierarchical_phenotype
    [pheno,parameters_list,pheno_non_zscore,pheno_with_outliers] = get_pheno_hierarchical(4, '/tmp/tmp.txt');
else
    [pheno,parameters_list,pheno_non_zscore,pheno_with_outliers] = get_pheno(4, '/tmp/tmp.txt');
end
if use_pyddm
    [pheno_not_hierarchical,parameters_list_not_hierarchical,~,~] = get_pheno(4, '/tmp/tmp.txt');
    idx = find(cellfun(@(x) contains(x,'rdm'), parameters_list_not_hierarchical));
    parameters_list(idx) = parameters_list_not_hierarchical(idx);
    pheno(idx,:,:) = pheno_not_hierarchical(idx,1:sub_num,:); % This works because it's the same number of parameters in both models
end


%% Calculate correlation between average/median value for each parameter (across weeks) and the different personality factors
corr_mat_pearson = nan(length(parameters_list),length(factor_names));
corr_mat_spearman = nan(length(parameters_list),length(factor_names));

for pI = 1:length(parameters_list)
    subjects_to_remove = [];
    tmp1 = nan(90,1);
    if strcmp(correlate_with_phenotype_mean_or_std,'mean')
        tmp1 = nanmean(squeeze(pheno(pI,:,:)),2);
    elseif strcmp(correlate_with_phenotype_mean_or_std,'median')
        tmp1 = nanmedian(squeeze(pheno(pI,:,:)),2);
    elseif strcmp(correlate_with_phenotype_mean_or_std,'std')
        tmp1 = nanstd(squeeze(pheno(pI,:,:)),[],2);
    elseif strcmp(correlate_with_phenotype_mean_or_std,'TotalDifference')
        % Get the difference between the last and first session (even if they are not weeks 12 and 1)
        tmp = squeeze(pheno(pI,:,:));
        for sI = 1:size(tmp,1)
            idx = find(~isnan(tmp(sI,:)));
            if isempty(idx) % This could happen if this is an excluded subject, like in gng
                subjects_to_remove(end+1) = sI;
                continue
            end
            tmp1(sI) = tmp(sI,idx(end))-tmp(sI,idx(1));
        end
    elseif strcmp(correlate_with_phenotype_mean_or_std,'baseline')
        % Get the difference between the last and first session (even if they are not weeks 12 and 1)
        tmp = squeeze(pheno(pI,:,:));
        for sI = 1:size(tmp,1)
            idx = find(~isnan(tmp(sI,:)));
            if isempty(idx) % This could happen if this is an excluded subject, like in GNG or LT
                subjects_to_remove(end+1) = sI;
                continue
            end
            tmp1(sI) = tmp(sI,idx(1));
        end
    end
    if ~isempty(subjects_to_remove)
        tmp1(subjects_to_remove,:) = [];
    end
    for fI = 1:length(factor_names)
        tmp2 = squeeze(data_factor(1:sub_num,fI));
        if ~isempty(subjects_to_remove)
            tmp2(subjects_to_remove,:) = [];
        end
        r_pearson = corr([tmp1(:),tmp2(:)],'Type','Pearson');
        r_pearson = r_pearson(2);
        corr_mat_pearson(pI,fI) = r_pearson;


        idx = ~isnan(tmp1) & ~isnan(tmp2);
        if ~any(idx)
            continue
        end
        r_spearman = corr([tmp1(idx),tmp2(idx)],'Type','Spearman');
        r_spearman = r_spearman(2);
        corr_mat_spearman(pI,fI) = r_spearman;
    end
end

%% Repeat the correlation calculation for random permutations of the subjects
tic
perm_num = 1000;
corr_mat_pearson_perm = nan(length(parameters_list),factor_num,perm_num);
corr_mat_spearman_perm = nan(length(parameters_list),factor_num,perm_num);

for permI = 1:perm_num
    perm = randperm(sub_num);

    for pI = 1:length(parameters_list)
        subjects_to_remove = [];

        if strcmp(correlate_with_phenotype_mean_or_std,'mean')
            tmp1 = nanmean(squeeze(pheno(pI,:,:)),2);
        elseif strcmp(correlate_with_phenotype_mean_or_std,'median')
            tmp1 = nanmedian(squeeze(pheno(pI,:,:)),2);
        elseif strcmp(correlate_with_phenotype_mean_or_std,'std')
            tmp1 = nanstd(squeeze(pheno(pI,:,:)),[],2);

        elseif strcmp(correlate_with_phenotype_mean_or_std,'TotalDifference')
            % Get the difference between the last and first session (even if they are not weeks 12 and 1)
            tmp = squeeze(pheno(pI,:,:));
            for sI = 1:size(tmp,1)
                idx = find(~isnan(tmp(sI,:)));
                if isempty(idx) % This could happen if this is an excluded subject, like in gng
                    subjects_to_remove(end+1) = sI;
                    continue
                end
                tmp1(sI) = tmp(sI,idx(end))-tmp(sI,idx(1));
            end
        elseif strcmp(correlate_with_phenotype_mean_or_std,'baseline')
        % Get the difference between the last and first session (even if they are not weeks 12 and 1)
        tmp = squeeze(pheno(pI,:,:));
        for sI = 1:size(tmp,1)
            idx = find(~isnan(tmp(sI,:)));
            if isempty(idx) % This could happen if this is an excluded subject, like in gng
                subjects_to_remove(end+1) = sI;
                continue
            end
            tmp1(sI) = tmp(sI,idx(1));
        end
    end
        if ~isempty(subjects_to_remove)
            tmp1(subjects_to_remove,:) = [];
        end

        for fI = 1:length(factor_names)
            tmp2 = squeeze(data_factor(1:sub_num,fI));
            if ~isempty(subjects_to_remove) 
                tmp2 = tmp2(perm(~ismember(perm,subjects_to_remove))); % Permute one of the variables (the factors) but not the other (the phenotype parameters)
            end
            
            r_pearson = corr([tmp1(:),tmp2(:)],'Type','Pearson');
            r_pearson = r_pearson(2);
            corr_mat_pearson_perm(pI,fI,permI) = r_pearson;

            idx = ~isnan(tmp1) & ~isnan(tmp2);
            if ~any(idx)
                continue
            end
            r_spearman = corr([tmp1(idx),tmp2(idx)],'Type','Spearman');
            r_spearman = r_spearman(2);
            corr_mat_spearman_perm(pI,fI,permI) = r_spearman;
        end
    end
end
toc


%% Calculate a different threshold per pair
multiple_comparisons_scaling_factor = 1;
for pI = 1:length(parameters_list)
    for fI  = 1:factor_num
        tmp = abs(corr_mat_pearson_perm(pI,fI,:));
        tmp = tmp(tmp>0);
        thresh_per_pair(pI,fI) = prctile(tmp,100-(alpha_statistic/multiple_comparisons_scaling_factor)*100);

        tmp = abs(corr_mat_spearman_perm(pI,fI,:));
        tmp = tmp(tmp>0);
        thresh_per_pair_spearman(pI,fI) = prctile(tmp,100-(alpha_statistic/multiple_comparisons_scaling_factor)*100);
    end
end
% Calculate threshold with correction for multiple comparisons
multiple_comparisons_scaling_factor = length(parameters_list)*factor_num; % Multiple comparisons affect the percentile of the data we want to exceed
for pI = 1:length(parameters_list)
    for fI  = 1:factor_num
        tmp = abs(corr_mat_pearson_perm(pI,fI,:));
        tmp = tmp(tmp>0);
        thresh_per_pair_multi_comp_corrected(pI,fI) = prctile(tmp,100-(alpha_statistic/multiple_comparisons_scaling_factor)*100);

        tmp = abs(corr_mat_spearman_perm(pI,fI,:));
        tmp = tmp(tmp>0);
        thresh_per_pair_spearman_multi_comp_corrected(pI,fI) = prctile(tmp,100-(alpha_statistic/multiple_comparisons_scaling_factor)*100);
    end
end


% Plot the correlations
figure('color','w')
subplot(1,2,1)
imagesc(corr_mat_pearson);
set(gca,'xtick',1:length(factor_names),'xticklabel',factor_names,'XTickLabelRotation',-30)
set(gca,'ytick',1:length(parameters_list),'yticklabel',strrep(parameters_list,'_',' '));
caxis([-1 1])
colormap(coolwarm)
title('Pearson')
if strcmp(correlate_with_phenotype_mean_or_std,'std')
    title('Pearson (corr with phenotype std)')
end
if strcmp(correlate_with_phenotype_mean_or_std,'median')
    title('Pearson (corr with phenotype median)')
end
if strcmp(correlate_with_phenotype_mean_or_std,'TotalDifference')
    title('Pearson (corr with phenotype change (week 12 - week 1)')
end
if strcmp(correlate_with_phenotype_mean_or_std,'baseline')
    title('Pearson (corr with phenotype baseline (week 1)')
end

set(gca,'fontsize',14)


for pI = 1:length(parameters_list)
    for fI = 1:factor_num
        % Add * for significant entries in the upper right half (Pearson
        % correlations)
        if abs(corr_mat_pearson(pI,fI))>(thresh_per_pair(pI,fI))
            text(fI,pI,'x','FontSize',14)
        end
        if abs(corr_mat_pearson(pI,fI))>(thresh_per_pair_multi_comp_corrected(pI,fI))
            text(fI,pI,'*','FontSize',14)
        end
    end
end

subplot(1,2,2)
imagesc(corr_mat_spearman);
set(gca,'xtick',1:length(factor_names),'xticklabel',factor_names,'XTickLabelRotation',-30)
set(gca,'ytick',1:length(parameters_list),'yticklabel',strrep(parameters_list,'_',' '));
caxis([-1 1])
colormap(coolwarm)
title('Spearman')
if strcmp(correlate_with_phenotype_mean_or_std,'std')
    title('Spearman (corr with phenotype std)')
end
if strcmp(correlate_with_phenotype_mean_or_std,'median')
    title('Spearman (corr with phenotype median)')
end
if strcmp(correlate_with_phenotype_mean_or_std,'TotalDifference')
    title('Spearman  (corr with phenotype change (week 12 - week 1)')
end
if strcmp(correlate_with_phenotype_mean_or_std,'baseline')
    title('Pearson (corr with phenotype baseline (week 1)')
end


set(gca,'fontsize',14)


for pI = 1:length(parameters_list)
    for fI = 1:factor_num
        % Add * for significant entries in the lower left half (Spearman
        % correlations)
        if abs(corr_mat_pearson(pI,fI))>(thresh_per_pair_spearman(pI,fI))
            text(fI-0.1,pI,'x','FontSize',14)
        end
        if abs(corr_mat_pearson(pI,fI))>(thresh_per_pair_spearman_multi_comp_corrected(pI,fI))
            text(fI,pI,'*','FontSize',14)
        end

    end
end

set(gcf,'position',get(0,'screensize'))

fig_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs';
if export_figures
    fig_file = fullfile(fig_dir,['hierarchical_completely_independent_parameter_personality_correlations_',num2str(factor_num),'_personalityFactors.png']);
    export_fig(gcf,fig_file,'-dpng','-r200');
end





%% A piece of code from a previous version where I explicitly extracted pi_gng (and didn't call it "param")
% pi_gng_sum = sum(pi_gng,2);
% sub_high_pi = find(pi_gng_sum>15);
%
% grouping =  zeros(size(subjects)); % Grouping for boxplot
% grouping(sub_high_pi) = 1;
%
% figure('color','w')
% for fI = 1:length(factor_names)
%     subplot(1,6,fI)
%     boxplot(data_factor(:,fI),grouping);
%
%
%     [h,p,ci,stats] = ttest2(data_factor(grouping==0,fI),data_factor(grouping==1,fI));
%
%     title({factor_names{fI}; ['p=',num2str(p)]})
%     axis square
% end
%
%
% figure('color','w')
% for fI = 1:length(factor_names)
%     subplot(1,6,fI)
%     h = histogram(data_factor(grouping==0,fI));
%     hold all
%     h = histogram(data_factor(grouping==1,fI),'binedges',h.BinEdges);
%     axis square
% end

%%
% for fI = 1:3
% [h,p,ci,stats] = ttest2(mean(exo_pca(fI,grouping==0,:),3),mean(exo_pca(fI,grouping==1,:),3))
% % if p<0.1
% %    disp(fI)
% % end
%
% end