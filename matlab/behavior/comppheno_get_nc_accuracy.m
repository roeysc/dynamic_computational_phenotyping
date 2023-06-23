function [accuracy_per_week,rt_per_week] = comppheno_get_nc_accuracy()
% Get accuracy for each subject in each one of the conditions, and in each
% session (subject x week x condition).

% Get subjects in the order we always use them
subjects = comppheno_get_subjects();

% Load the CD data
fname = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/nc_data_for_stan_90s.csv';
data = readtable(fname);

%%
accuracy_per_week= [];
for sI = 1:length(subjects)
    sub = subjects{sI};
    for wI = 1:12
        idx = (strcmp(data.subjectId,sub) & data.weekId==wI);
        if all(idx==0)
            accuracy_per_week(sI,wI) = nan;
            rt_per_week(sI,wI) = nan;
            continue
        end
        idx_correct = (strcmp(data.subjectId,sub) & data.weekId==wI) & strcmp(data.correct,'True');
        accuracy_per_week(sI,wI) = sum(idx_correct)/sum(idx);
        % Get RT
        rt_per_week(sI,wI) = nanmean(data.rt(idx));

    end
end