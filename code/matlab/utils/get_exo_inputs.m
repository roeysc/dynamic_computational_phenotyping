function exo = get_exo_inputs(task)
% Returns the exogenous variables matrix (35 PCs x 88 subjects x 12 weeks).
% This version of the function also allows to calculate
%
% Inputs:
% task: string of the task ('lottery','itc', etc.). If given, we will use
%   the specific life survey data collected on the day of this task, rather
%   than using the average data across all 3 days in this week.

% Parse input arguments.
if nargin==0
    task = '';
end


comppheno_set_dirs % Load the comppheno_dir variable

    'cd_data_for_stan_90s.csv');

if strcmp(task,'gng') || strcmp(task,'rdm')
    exo_file = fullfile(comppheno_data_dir,'survey_pca','exo_big_data_day1_gng_rdm.csv';
elseif strcmp(task,'nc') || strcmp(task,'smb') || strcmp(task,'lottery')
    exo_file = fullfile(comppheno_data_dir,'survey_pca','exo_big_data_day2_nc_smb_lt.csv';
elseif strcmp(task,'cd') || strcmp(task,'itc')
    exo_file = fullfile(comppheno_data_dir,'survey_pca','exo_big_data_day3_cd_itc.csv';
else
    exo_file = fullfile(comppheno_data_dir,'survey_pca','exo_data_for_stan.csv';
end

subjects = comp_pheno_get_subjects(); % Do not use the input argument 'csv_file' here, as some exo files have more than 90 subjects, and we just want our 90 subjects in their correct order
data = readtable(exo_file);
tmp = regexp(data.Properties.VariableNames,'^q'); % column indices matching questions in the table
qidx = find(~cellfun(@isempty,tmp));
qnum = length(qidx); % Number of question (exogenous inputs)

for sI = 1:length(subjects)
   sub = subjects{sI};
       for wI = 1:12
          for qI = 1:qnum
           idx = find(strcmp(data.subjectId,sub) & data.weekId==wI);
           if ~isempty(idx)
               exo(qI,sI,wI) = table2array(data(idx(1),qidx(qI))); % idx(1) is used just because one subject has a duplicate week (subject 90, week 11)
           else
               exo(qI,sI,wI) = nan;
           end
       end
   end
end


% Fill in missing sessions with the average of existing sessions in the
% same week
if ~isempty(task)
    exo_average = get_exo_inputs(); % Giving no input arguments means getting the average exogenous variables from all 3 days of each week (or less, if some sessions were missing)
    idx = find(isnan(exo));
    exo(idx) = exo_average(idx);
end

%% z-score each input variable
% for qI = 1:size(exo,1)
% x = exo(qI,:,:);
% x = x(:);
% x_zscore = (x - nanmean(x))/nanstd(x);
% tmp = reshape(x_zscore,size(exo,2),size(exo,3));
% exo(qI,:,:) = tmp;
% end

%% Check plots
% figure('color','w')
% subplot(1,3,1)
% imagesc(squeeze(exo1(1,:,:)));
% subplot(1,3,2)
% imagesc(squeeze(exo2(1,:,:)));
% subplot(1,3,3)
% imagesc(squeeze(exo3(1,:,:)));

% % Scatter plots of the exo variables across days, to visualize the
% % variability between days.
% for qI = 1:35
%     tmp1 = squeeze(exo1(qI,:,:));
%     tmp2 = squeeze(exo2(qI,:,:));
%     tmp3 = squeeze(exo3(qI,:,:));
%     figure('color','w')
%     
%     scatter(tmp1(:),tmp2(:),'filled','MarkerFaceAlpha',0.3);
%     hold all
%     scatter(tmp1(:),tmp3(:),'filled','MarkerFaceAlpha',0.3);
%     identityLine(gca)
%     title(num2str(qI))
% end
