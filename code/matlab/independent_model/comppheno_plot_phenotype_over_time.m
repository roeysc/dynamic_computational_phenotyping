addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/roey_code/dynamical_system/mle')
addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/roey_code/pareto')
addpath('/n/home10/rschurr/Matlab_code/external_toolboxes/altmany-export_fig-3e1e29a'); % For exporting figures


addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/publish_code/matlab')
clear

%% Settings
matplotlib_colormaps; % Load colormaps
use_hierarchical_phenotype = true;
use_pyddm = false; % Replace the DDM from Stan with pyddm
add_random_subjects = false;

correlate_with_phenotype_mean_or_std = 'median'; % mean, median or std
factor_num = 3; % Number of personality factors to use (3 or 6)
sub_num = 90;
model = 'independent'; % ***independent / super_model

use_exo_per_task = false;

matplotlib_colormaps
add_power_law_fit = false;
export_figures = true;

%%
for use_pheno_zscored = 0%0:1 % *** Notice this
    if use_pheno_zscored
        fig_dir = fullfile('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs/phenotype_over_time_zscored',model);
    else
        fig_dir = fullfile('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs/phenotype_over_time',model);
        fig_dir = fullfile('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs/phenotype_over_time/unmodified_gng_model/',model);
    end
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir)
    end
if ~exist(fullfile(fig_dir,'pdf'),'dir')
    mkdir(fullfile(fig_dir,'pdf'));
end

    %% Get the phenotype
    if use_hierarchical_phenotype
        [pheno,parameters_list,pheno_non_zscore,pheno_with_outliers] = get_pheno_hierarchical(4, '/tmp/tmp2.txt',model);
    else
        [pheno,parameters_list,pheno_non_zscore,pheno_with_outliers] = get_pheno(4, '/tmp/tmp.txt');
    end
    if use_pyddm
        [pheno_not_hierarchical,parameters_list_not_hierarchical,~,~] = get_pheno(4, '/tmp/tmp.txt');
        idx = find(cellfun(@(x) contains(x,'rdm'), parameters_list_not_hierarchical));
        parameters_list(idx) = parameters_list_not_hierarchical(idx);
        pheno(idx,:,:) = pheno_not_hierarchical(idx,1:sub_num,:); % This works because it's the same number of parameters in both models
    end

    %% Get the exogenous variables
    max_pc = 3; % Take up to 3 PCs from the exogenous vairables
    if use_exo_per_task
        session_to_hold_out = []; % Don't hold out any sessions
        data_exo_pca_day1 = get_exo_inputs_after_PCA(session_to_hold_out,'gng');
        data_exo_pca_day1 = data_exo_pca_day1(1:max_pc,:,:);
        data_exo_pca_day2 = get_exo_inputs_after_PCA(session_to_hold_out,'nc');
        data_exo_pca_day2 = data_exo_pca_day2(1:max_pc,:,:);
        data_exo_pca_day3 = get_exo_inputs_after_PCA(session_to_hold_out,'cd');
        data_exo_pca_day3 = data_exo_pca_day3(1:max_pc,:,:);
    else
        task = ''; % Use the average across days to have all tasks. &&&TODO: Actually, it might be better to use the relevant exo for each task
        session_to_hold_out = []; % Don't hold out any sessions
        data_exo_pca = get_exo_inputs_after_PCA(session_to_hold_out,task); % This will exclude the current session from the PCA calculation. It will also use the 'task' variable to know which day to take the exogenous variables from (day 1/2/3 of the week, or their average if 'task' is empty).
        data_exo_pca = data_exo_pca(1:max_pc,:,:);
    end

    % Combine the exogenous variables with the phenotype
    pheno_non_zscore = cat(1,pheno_non_zscore,data_exo_pca);
    parameters_list{end+1} = 'Valence';
    parameters_list{end+1} = 'Arousal';
    parameters_list{end+1} = 'Self control';

    %% Plot the phoenotype over time
    for pI = 1:length(parameters_list)

        if use_pheno_zscored
            param = squeeze(pheno(pI,:,:));
        else
            param = squeeze(pheno_non_zscore(pI,:,:));
        end
        param_name = strrep(parameters_list{pI},'_',' ');

        exclude_subject_choosing_safe_more_than_80_percent = true;
        if exclude_subject_choosing_safe_more_than_80_percent && any(strfind(parameters_list{pI},'_lt'))
            idx = [1,4,7,8,9,11,14,15,19,23,24,25,28,32,35,36,37,47,50,53,58,60,61,65,66,68,70,74,77,83,84,88];
            param(idx,:) = [];
        end
        %% Parameter heat map
        tmp = param;
        %         figure('color','w')
        %         nanimagesc(tmp)
        %         caxis(prctile(tmp(:),[2,98]))
        %         colorbar
        %         set(gca,'tickdir','out','fontsize',14,'box','off','xtick',1:12)
        %         xlabel('session')
        %         ylabel('subject')
        %         title(param_name)
        %         axis square
        %
        %         set(gcf,'position',get(0,'screensize'))
        %         if export_figures
        %             fig_file = fullfile(fig_dir,sprintf('hierarchical_%s_%s_heatmap.png',model,param_name));
        %             export_fig(gcf,fig_file,'-dpng','-r200');
        %         end

        %% Parameter curve
        figure('color','w')
        %         h = plot(nanmean(tmp),'.k','linew',2);
        h = scatter(1:12,nanmean(tmp),128,[0,0,0],'filled');
        hold all
        herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color','k','LineStyle','none','LineWidth',2)
        ylims = ylim;
        ylims(1) = floor(ylims(1)*10)/10;
        ylims(2) = ceil(ylims(2)*10)/10;
        ylim(ylims); % Round accuracy values to closest 10,20, etc.
        set(gca,'FontSize',24,'TickDir','out','linewidth',2,'box','off','xtick',1:12)
        xlabel('session')
        ylabel(param_name)
        title(param_name)
        axis square
        % h = plot(tmp','linew',0.5,'Color',0.5*[1,1,1]); % Add individual subjects
        pause(1)
        if add_random_subjects
            hfig = gcf;
            idx = randperm(size(tmp,1));
            hfigtmp = figure;
            hold all
            for sI = 1:10
                h = plot(tmp(idx(sI),:)-nanmean(tmp(idx(sI),:)),'-','linew',2);
            end
            ylims = ylim;
            close(gcf);
            figure(hfig);
            for sI = 1:10
                h = plot(tmp(idx(sI),:)-nanmean(tmp(idx(sI),:)),'-','linew',2);
            end
            ylim(ylims);
        end
        %%
        if add_power_law_fit
            x = 1:12;
            y = nanmean(tmp);
            [xData, yData] = prepareCurveData( x, y );
            % Set up fittype and options.
            ft = fittype( 'power1' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display='off';
            opts.StartPoint = [0.455781095951251 0.170138339922848];
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            %hfit = plot(fitresult);
            t = 1:0.1:12;
            hfit = plot(t,fitresult.a.*t.^fitresult.b);
            hfit.Color = h.Color;%0.5*[1,1,1]
            hfit.LineWidth = 2;
            subtitle(['R^2_a_d_j = ', num2str(gof.adjrsquare,2)])
        end
        set(gcf,'position',get(0,'screensize'))

        if export_figures
            set(gcf,'position',get(0,'ScreenSize'));
            fig_file = fullfile(fig_dir,sprintf('hierarchical_%s_%s_over_time.png',model,param_name));
            export_fig(gcf,fig_file,'-dpng','-r200');
            fig_file = fullfile(fig_dir,'pdf',sprintf('hierarchical_%s_%s_over_time.pdf',model,param_name));
            exportgraphics(gcf,fig_file,'ContentType','vector')
        end

        pause(0.5)
        %         clear h
        %         clear hfit

        %         close all
    end




end