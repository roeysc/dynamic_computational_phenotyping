function [accuracy_per_week_and_coherence,rt_per_week_and_coherence] = comppheno_get_rdm_accuracy()
% Get accuracy for each subject in each one of the conditions, and in each
% session (subject x week x coherence level).

% Get subjects in the order we always use them
subjects = comppheno_get_subjects();

% Load the RDM data
comppheno_set_dirs % Load the comppheno_dir variable
fname = fullfile(comppheno_data_dir,'rdm_data_for_stan_90s.csv');
data = readtable(fname);

accuracy = [];
coherence_vals = [0.05,0.1,0.35,0.5];
for sI = 1:length(subjects)
    sub = subjects{sI};
    for cond = 1:length(coherence_vals) % 4 coherence values
        coh = coherence_vals(cond);
        for wI = 1:12
            idx = (strcmp(data.subjectId,sub) & data.coh==coh & data.weekId==wI);
            if all(idx==0)
                accuracy_per_week_and_coherence(sI,wI,cond) = nan;
                rt_per_week_and_coherence(sI,wI,cond) = nan;
                continue
            end
            idx_correct = (strcmp(data.subjectId,sub) & data.coh==coh& data.weekId==wI & strcmp(data.correct,'True'));
            accuracy_per_week_and_coherence(sI,wI,cond) = sum(idx_correct)/sum(idx);
            % Get RT
            rt_per_week_and_coherence(sI,wI,cond) = nanmean(data.rt(idx));
        end 
    end
end