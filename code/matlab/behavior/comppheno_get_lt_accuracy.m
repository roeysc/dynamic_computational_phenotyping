function [accuracy_per_week_and_cond,rt_per_week_and_cond] = comppheno_get_lt_accuracy()
% Get accuracy for each subject in each one of the conditions, and in each
% session (subject x week x cond level).

% Get subjects in the order we always use them
subjects = comppheno_get_subjects();

% Load the lt data
comppheno_set_dirs % Load the comppheno_dir variable
fname = fullfile(comppheno_data_dir,'lt_data_for_stan_90s.csv');
data = readtable(fname);

accuracy = [];
conditions = {'small','medium','large'};
for sI = 1:length(subjects)
    sub = subjects{sI};

for wI = 1:12
    idx = (strcmp(data.subjectId,sub) & data.weekId==wI);
        if all(idx==0)
            accuracy_per_week_and_cond(sI,wI,cond) = nan;
            rt_per_week_and_cond(sI,wI,cond) = nan;
            continue
        end
    
    % Low rewards    
    cond = 1;
    idx = (strcmp(data.subjectId,sub) & data.hi_narr<10 & data.weekId==wI);
    idx_correct = strcmp(data.subjectId,sub) & data.hi_narr<10 & data.weekId==wI & strcmpi(data.correct,'True');
    accuracy_per_week_and_cond(sI,wI,cond) = sum(idx_correct)/sum(idx);
    rt_per_week_and_cond(sI,wI,cond) = nanmean(data.rt(idx));

    % Medium rewards
    cond = 2;
    idx = (strcmp(data.subjectId,sub) & data.hi_narr>10  & data.hi_narr<100 & data.weekId==wI);
    idx_correct = strcmp(data.subjectId,sub) & data.hi_narr>10  & data.hi_narr<100 & data.weekId==wI & strcmpi(data.correct,'True');
    accuracy_per_week_and_cond(sI,wI,cond) = sum(idx_correct)/sum(idx);
    rt_per_week_and_cond(sI,wI,cond) = nanmean(data.rt(idx));

    % High rewards    
    cond = 3;
    idx = (strcmp(data.subjectId,sub) & data.hi_narr>100 & data.weekId==wI);
    idx_correct = strcmp(data.subjectId,sub) & data.hi_narr>100 & data.weekId==wI & strcmpi(data.correct,'True');
    accuracy_per_week_and_cond(sI,wI,cond) = sum(idx_correct)/sum(idx);
    rt_per_week_and_cond(sI,wI,cond) = nanmean(data.rt(idx));



end
end

%%
choice_std=nan(90,12);

for sI = 1:length(subjects)
    sub = subjects{sI};

for wI = 1:12
    idx = (strcmp(data.subjectId,sub) & data.weekId==wI);
        if all(idx==0)
            continue
        end
choice_std(sI,wI)= std(data.choice(idx));    
end
end
