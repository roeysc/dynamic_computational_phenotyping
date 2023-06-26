function [rt_per_week] = comppheno_get_itc_accuracy()
% Get RT for each subject and session (subject x week).
% Since we don't have a way of estimating accuracy in this task, we only get RT here (the naming of the funciton remained"accuracy" for consistency with other tasks).

% Get subjects in the order we always use them
subjects = comppheno_get_subjects();

% Load the itc data
comppheno_set_dirs % Load the comppheno_dir variable
fname = fullfile(comppheno_data_dir,'itc_data_for_stan_90s.csv');
data = readtable(fname);

rt_per_week = [];
for sI = 1:length(subjects)
    sub = subjects{sI};

for wI = 1:12
    idx = (strcmp(data.subjectId,sub) & data.weekId==wI);
        if all(idx==0)
            rt_per_week(sI,wI) = nan;
            continue
        end
    
    % Low rewards    
    rt_per_week(sI,wI) = nanmean(data.rt(idx));

end
end