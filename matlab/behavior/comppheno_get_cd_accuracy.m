function [accuracy_per_week_and_condition, rt_per_week_and_condition,rt_per_week_and_condition_prctile,rt_per_week_and_condition_max,rt_per_week_and_condition_outlier] = comppheno_get_cd_accuracy()
% Get accuracy for each subject in each one of the conditions, and in each
% session (subject x week x condition).

% Get subjects in the order we always use them
subjects = comppheno_get_subjects();

% Load the CD data
fname = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/cd_data_for_stan_90s.csv';
fname = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/cd_data_for_stan_90s_with_rt.csv';

data = readtable(fname);

blocks = {'three_t1','four_t1','six_t1','eight_t1','eight_tvary'};

accuracy_per_week_and_condition = [];
rt_per_week_and_condition = [];
rt_per_week_and_condition_prctile = [];
rt_per_week_and_condition_max = [];
rt_per_week_and_condition_outlier = [];

for sI = 1:length(subjects)
    sub = subjects{sI};
    for wI = 1:12
        for bI = 1:length(blocks)
            idx = (strcmp(data.subjectId,sub) & strcmp(data.block_condition,blocks{bI}) & data.weekId==wI & data.rt<10*1000);  % Ignore trials with RT>10 seconds

            if isempty(idx) | sum(idx)==0
                accuracy_per_week_and_condition(sI,wI,bI) = nan;
                rt_per_week_and_condition(sI,wI,bI) = nan;
                rt_per_week_and_condition_max(sI,wI,bI) = nan;
            rt_per_week_and_condition_outlier(sI,wI,bI) = nan;
                continue
            end
            % Get RT
            rt_per_week_and_condition(sI,wI,bI) = nanmean(data.rt(idx));
            rt_per_week_and_condition_prctile(sI,wI,bI) = prctile(data.rt(idx),90);
            rt_per_week_and_condition_max(sI,wI,bI) = max(data.rt(idx));
            rt_per_week_and_condition_outlier(sI,wI,bI) = sum(data.rt(idx)>10*1000); % Max 10 seconds
    
%             continue % TODO: Fix this when we have the corrected final CD data.
            idx_correct = (strcmp(data.subjectId,sub) & strcmp(data.block_condition,blocks{bI}) & data.weekId==wI)  & data.rt<10*1000 & (logical(data.t)==logical(data.behav)); % The last condition makes sure that even cases with four targets (data.t=4) will be considered correct given a behavior of data.behav=1 (subject indicated a change)
            accuracy_per_week_and_condition(sI,wI,bI) = sum(idx_correct)/sum(idx);

        end
        
    end
end