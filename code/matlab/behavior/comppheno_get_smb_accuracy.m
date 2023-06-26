function [accuracy_per_week_and_condition,rt_per_week_and_condition,reward_per_week_and_condition,accuracy_per_week_and_condition_per_half] = comppheno_get_smb_accuracy()
% Get reward for each subject in each one of the conditions, and in each
% session (subject x week x condition).

% Get subjects in the order we always use them
subjects = comppheno_get_subjects();


%%
% Load the SMB data
comppheno_set_dirs % Load the comppheno_dir variable
fname = fullfile(comppheno_data_dir,'smb_data_for_stan_90s.csv');
data = readtable(fname);

reward_per_week_and_condition = [];
rt_per_week_and_condition = [];
accuracy_per_week_and_condition = []; % Accuracy defined here as choosing the arm with the higher mean: RS, SR, RR, SS
accuracy_per_week_and_condition_per_half = []; % Accuracy also for each third (first 5 trials, last 5 trials)
for sI = 1:length(subjects)
    sub = subjects{sI};
    for wI = 1:12
        for cond = 1:4
            idx = (strcmp(data.subjectId,sub) & data.cond==cond & data.weekId==wI);
            if isempty(idx)
                reward_per_week_and_condition(sI,wI,cond) = nan;
                accuracy_per_week_and_condition(sI,wI,cond) = nan;
                rt_per_week_and_condition(sI,wI,cond) = nan;
                accuracy_per_week_and_condition_per_half(sI,wI,cond) = nan;
                continue
            end
            reward_per_week_and_condition(sI,wI,cond) = sum(data.reward(idx));
            rt_per_week_and_condition(sI,wI,cond) = nanmean(data.rt(idx));
        
            tmp = []; % hold the accuracy of each block, to average them later
            tmp2 = []; % hold the accuracy of each block (only second half), to average them later
            for bI = 1:30
                idx = find(strcmp(data.subjectId,sub) & data.cond==cond & data.weekId==wI & strcmp(data.block,['block_',num2str(bI)]));
                if isempty(idx) || sum(idx)==0
                    continue
                end
                m1 = data.mean1(idx(1));
                m2 = data.mean2(idx(1));
                if m1>m2
                    tmp(end+1) = sum(data.chosen_machine(idx)==1)./length(idx); % Best machine (1) was chosen
                    tmp2(end+1) = sum(data.chosen_machine(idx(6:end))==1)./length(idx(6:end)); % Best machine (1) was chosen
                else
                    tmp(end+1) = sum(data.chosen_machine(idx)==0)./length(idx); % Best machine (2) was chosen
                    tmp2(end+1) = sum(data.chosen_machine(idx(6:end))==0)./length(idx(6:end)); % Best machine (2) was chosen
                end
            end
             accuracy_per_week_and_condition(sI,wI,cond) = mean(tmp);
             accuracy_per_week_and_condition_per_half(sI,wI,cond) = mean(tmp2);
        end
        
    end
end