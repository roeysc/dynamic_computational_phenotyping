function [accuracy,accuracy_per_week, accuracy_per_week_and_block,accuracy_per_week_and_block_second_half,rt_per_week] = comppheno_get_gng_accuracy()
% Get accuracy for each subject in each one of the conditions, and in each
% session.

% Get subjects in the order we always use them
subjects = comppheno_get_subjects();

% Load the GNG data
comppheno_set_dirs % Load the comppheno_dir variable
fname = fullfile(comppheno_data_dir,'gng_data_for_stan_90s.csv');

data = readtable(fname);

accuracy = [];
for sI = 1:length(subjects)
    sub = subjects{sI};
    for cond = 1:4 % conditions = {'Go2win','NoGo2win ','Go2avoid','NoGo2avoid'};
        % Accuracy across all sessions
        idx = (strcmp(data.subjectId,sub) & data.cond==cond);
        idx_correct = (strcmp(data.subjectId,sub) & data.cond==cond & strcmp(data.correct,'True'));
        accuracy(sI,cond) = sum(idx_correct)/sum(idx);
        
        % Accuracy per session
        for wI = 1:12
            idx = (strcmp(data.subjectId,sub) & data.cond==cond & data.weekId==wI);
            if all(idx==0)
                accuracy_per_week(sI,wI,cond) = nan;
                rt_per_week(sI,wI,cond) = nan;
                continue
            end
            idx_correct = (strcmp(data.subjectId,sub) & data.cond==cond & data.weekId==wI & strcmp(data.correct,'True'));
            accuracy_per_week(sI,wI,cond) = sum(idx_correct)/sum(idx);
            rt_per_week(sI,wI,cond) = nanmean(data.rt(idx));
        end
        
        % Accuracy per session and block (to see if maye divergent sessions
        % are sessions with a large discrepancy between different blocks)
        for wI = 1:12
            for bI = 1:3
                block_str = ['block_' num2str(bI)];
                idx = (strcmp(data.subjectId,sub) & data.cond==cond & data.weekId==wI & strcmp(data.block,block_str));
                if all(idx==0)
                    accuracy_per_week_and_block(sI,wI,cond,bI) = nan;
                    accuracy_per_week_and_block_second_half(sI,wI,cond,bI) = nan;
                    continue
                end
                idx_correct = (strcmp(data.subjectId,sub) & data.cond==cond & data.weekId==wI & strcmp(data.block,block_str) & strcmp(data.correct,'True'));
                accuracy_per_week_and_block(sI,wI,cond,bI) = sum(idx_correct)/sum(idx);
                % Second half of each block and condition
                idx_tmp = find(idx);
                idx_tmp = idx_tmp(ceil(length(idx_tmp)/2):end);
                accuracy_per_week_and_block_second_half(sI,wI,cond,bI) = sum(strcmp(data.correct(idx_tmp),'True'))/length(idx_tmp);
            end
            end
        
    end
end


%% Learners and non-learners, following the original paper (See Supplementary there)
% overall_learners = mean(accuracy_per_week,3)>0.6;
% accuracy_last_block = squeeze(accuracy_per_week_and_block(:,:,:,2));
% 
% end_of_experiment_learners = sum(accuracy_last_block>0.8,3)==4;
% 
% learners = overall_learners & end_of_experiment_learners;
% 
% figure('color','w')
% histogram(sum(learners'))
% xlabel('"learnt" sessions')
% ylabel('Subject count')
% box off
% axis square