function exo_pca = comppheno_get_exo_inputs_after_PCA(session_to_hold_out,task)
% Returns the exogenous variables matrix (35 PCs x 88 subjects x 12 weeks)
% after PCA. Based on the plots below, taking the first 3 components sounds
% like a good idea.
% This version of the function also allows to calculate
% the PCA while witholding a single session out. This is useful for
% cross-validation of the Gaussian-Process and other models of the
% phenotype dynamics. Holding out a session in the PCA makes sure that we
% don't have data from this sessions "leaking" into our analysis prior to
% testing it on our held-out session. (The only place where data does
% "leak" is when we z-score prior to calculating the PCA. But that's
% reasonable. If we don't include the held-out session in the z-scoring,
% it's hard to project it to the newly calculated PCA basis.
% The held-out data is then projected to the PCA basis to be included in the
% function output, exo_pca.
%
% Inputs:
% session_to_hold_out: 1x2 vector of subject-week pair to remove when
%   calculating the PCA. This is useful if we want to calculate PCA on a
%   all sessions excluding one (and then projecting this one to the new
%   axis defined by the PCA)
% task: string of the task ('lottery','itc', etc.). If given, we will use
%   the specific life survey data collected on the day of this task, rather
%   than using the average data across all 3 days in this week.
%
% For example, calculate PCA while witholding subejct 3, week 12, and using only
% the data collected on the same day as the Random Dot Motion task:
% get_exo_inputs_after_PCA([3,12],'rdm')

% Parse input arguments.
if nargin==0
    session_to_hold_out = [];
    task = '';
elseif nargin==1
    task = '';
end

% % % if strcmp(task,'gng') || strcmp(task,'rdm')
% % %     exo_file = '/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/exo_big_data_day1_gng_rdm.csv';
% % % elseif strcmp(task,'nc') || strcmp(task,'smb') || strcmp(task,'lottery')
% % %     exo_file = '/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/exo_big_data_day2_nc_smb_lt.csv';
% % % elseif strcmp(task,'cd') || strcmp(task,'itc')
% % %     exo_file = '/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/exo_big_data_day3_cd_itc.csv';
% % % else
% % %     exo_file = '/n/gershman_lab/users/Stan/2_task_data_files/exo_data_for_stan.csv';
% % % end
% % % exo = get_exo_inputs(exo_file);
% % % exo_orig = exo;
% % %
% % % % Fill in missing sessions with the average of existing sessions in the
% % % % same week
% % % if ~isempty(task)
% % %     exo_file_average = '/n/gershman_lab/users/Stan/2_task_data_files/exo_data_for_stan.csv';
% % %     exo_average = get_exo_inputs(exo_file_average);
% % %     idx = find(isnan(exo));
% % %     exo(idx) = exo_average(idx);
% % % end

% First get the exogenous variables from the regular function without PCA
exo = get_exo_inputs(task);

% % % drug_users = [3     5     7     9    13    15    19    21    22    23    27    35    38    40    41    43    44    45    46    47    48    49    56    58    59    62    66    68 69    70    71    72    74    77    78    81    82    83    85    89    90];
% % % non_drug_users = 1:90;
% % % non_drug_users (ismember(non_drug_users ,drug_users))=[];
% % % exo = exo(:,non_drug_users ,:); % *** Have only drug useres
% % % exo(32,:,:) = [];
% % % warning('ASASDASDASDASDASDASDAS')

size_exo_orig = size(exo);

%%%% *** Find drug users
% % % tmp = squeeze(exo(32,:,:));
% % % drugs = find(nansum(tmp')>0)
%%%% ***

% Expand to a question X sessions matrix
exo = exo(:,:);
% % % exo = exo([1:20,22:27],:);% *** Remove questions with binary answeres
% % % size_exo_orig = [26,90,12]; % *** Update the size
% % % warning('ASASDASDASDASDASDASDAS')
n_sessions = size(exo,2);
n_params = size(exo,1);

% Z-score across subjects
% Here we z-score explicitly, because we want to preserve the missing data
% here, so it's easier to hold out a single session below
% *TODO*: Is it better to z-score within subject, because we care about the individual's dynamics?
exo_mean = nanmean(exo,2);
exo_std = nanstd(exo,[],2);
exo = (exo-exo_mean)./exo_std;

% Remove missing sessions
idx_missing_columns = find(isnan(sum(exo)));
exo(:,idx_missing_columns) = [];


% Remove unwanted session, if any
if ~isempty(session_to_hold_out)
    tmp = zeros(size_exo_orig(2),size_exo_orig(3)); %  subjects x weeks
    tmp(session_to_hold_out(1),session_to_hold_out(2)) = 1;
    idx = find(tmp(:));
    exo_held_out = exo(:,idx);
    exo(:,idx) = [];
    idx_missing_columns(end+1) = idx; % Add this held out session to the list of missing sessions
end

% Run PCA
[coeff, score, latent, tsquared, explained, mu] = pca(exo');

% Change the sign of the first 3 PCs to consistently map valence, arousal and self-control
% Valence
if coeff(1,1)>0 % Feeling nervous should be negative, so that valence is positive for positive feelings
    coeff(:,1) = -coeff(:,1);
    score(:,1) = -score(:,1);
end
% Arousal
if coeff(6,2)<0 % Feeling determined should be positive and high
    coeff(:,2) = -coeff(:,2);
    score(:,2) = -score(:,2);
end
% Self-control
if coeff(33,3)>0 % Much use of drugs should be negative for high self control
     coeff(:,3) = -coeff(:,3);
     score(:,3) = -score(:,3);
end



% Impute missing sessions with NaNs
count = 1;
exo_pca = nan(n_sessions,n_params); % This has the size of the original n_sessions, including missing ones
% exo_not_pca = nan(n_sessions,n_params); % This uis used just to make sure we reshaped the exo_pca correctly

for session = 1:n_sessions
    if ismember(session,idx_missing_columns)
        continue
    end
    exo_pca(session,:) = score(count,:);
    %     exo_not_pca(session,:) = exo(:,count);
    count = count+1;
end

exo_pca = transpose(exo_pca);
exo_pca = reshape(exo_pca,size_exo_orig);


% Project held-out session to the PCA space (see
% https://www.mathworks.com/matlabcentral/answers/53259-how-to-project-a-new-point-to-pca-new-basis)
if ~isempty(session_to_hold_out)
    exo_held_out_score = exo_held_out'/coeff';
    exo_pca(:,session_to_hold_out(1),session_to_hold_out(2)) = exo_held_out_score';
end

% exo_not_pca = transpose(exo_not_pca);
% exo_not_pca = reshape(exo_not_pca,size(exo_orig));
% exo_not_pca = exo_not_pca(:,:);
% exo_not_pca(:,idx_missing_columns) = [];

% Make sure we reshaped everything correctly

%% Plots for sanity checks
plot_sanity_checks = true;
if plot_sanity_checks
    reorder_questions = true;
    % Import question names
    exo_file = '/n/gershman_lab/users/Stan/2_task_data_files/exo_data_for_stan.csv';

    t = readtable(exo_file);
    questions = t.Properties.VariableNames(6:end);
    questions = questions(1:35);
    for qI =1:length(questions)
        idx = strfind(questions{qI},'_');
        questions{qI} = questions{qI}(idx+1:end);
        questions{qI} = strrep(questions{qI},'_',' ');
    end

    if reorder_questions
        idx = [2,5,6,8,10,12,1,3,4,7,9,11,13:35];
        questions = questions(idx);
        coeff = coeff(idx,:);
    end
    figure('color','w')
    plot(explained,'-O')
    xlabel('PC')
    ylabel('% Var explained')

    figure('color','w')
    plot(coeff(:,1),'linew',3)
    hold all; plot(coeff(:,2),'linew',3)
    hold all; plot(coeff(:,3),'linew',3)
    set(gca,'xtick',1:length(questions))
    set(gca,'xticklabel',questions,'XTickLabelRotation',-30)
    ylabel('Loadings')
    box off
    legend({'PC1','PC2','PC3'})
    legend({'Valence','Arousal','Self control'})
    plot(get(gca,'xlim'),[0 0],'--k')
    set(gcf,'Position',get(0,'screensize'));
    set(gca,'fontsize',18,'tickdir','out','box','off')

    %% Plot with bars

    figure('color','w')
    for pI = 1:2%3
        subplot(3,1,pI)
        bar(coeff(:,pI),'linew',3)
        box off
        set(gca,'fontsize',18,'tickdir','out','box','off')
        set(gca,'xtick',1:length(questions))
        set(gca,'xticklabel','');
        if pI==2
            ylabel('Loadings')
        end
    end
    set(gca,'xticklabel',questions,'XTickLabelRotation',-30)
    set(gcf,'Position',get(0,'screensize'));
    %%
export_figures = true;

addpath('/n/home10/rschurr/Matlab_code/external_toolboxes/altmany-export_fig-3e1e29a'); % For exporting figures
fig_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs/exo_loadings/';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end
if export_figures
    fig_file = fullfile(fig_dir,['exo_loadings.png']);
    export_fig(gcf,fig_file,'-dpng','-r200');
     fig_file = strrep(fig_file,'.png','.pdf');
     exportgraphics(gcf,fig_file,'ContentType','vector')
end


end
