
addpath('/n/holyscratch01/gershman_lab/Lab/shared_code')
addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/shared_code/required_toolboxes/spiderplot');
addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/shared_code/super_model_dynamics_figures/')
addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/shared_code/required_toolboxes'); % For the HPDI.m function
addpath('/n/home10/rschurr/Matlab_code/external_toolboxes/altmany-export_fig-3e1e29a'); % For exporting figures
rmpath('/n/home10/rschurr/Matlab_code/external_toolboxes/export_fig-master/')
addpath('/n/home10/rschurr/Matlab_code/comp_pheno/')
% addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/shared_code/required_toolboxes/invprctile')
% addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/shared_code/required_toolboxes/boxplotGroup_vs_2.0.0')
addpath('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/shared_code/required_toolboxes')
%matplotlib_colormaps; % Load colormaps
x_prev = 0;
xticks_all = [];
params_all = [];

%%

in_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/4_hierarchical/R_files_for_plots/models_with_no_weights_and_effects_around_zero/exo2';

fig_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs/relative_STD_d';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
if ~exist(fullfile(fig_dir,'pdf'),'dir')
    mkdir(fullfile(fig_dir,'pdf'));
end
export_figures = 0;
 write_all_pd = 0;
 add_iqr = 0;
export_figures_whole_task = 0;
export_learning_figures = 0;
examine_exo_variability = false;
tasks = {'cd','nc','itc','lt','smb','gng','rdm'};

%% Get Exo to calculate exo variability across time
if examine_exo_variability
    % based on all existing data points
    max_pc = 3; % Take up to 3 PCs from the exogenous vairables
    use_exo_per_task = true;
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

    % Concatenate the data from different days of the week alongside the
    % weeks, to have all the datapoints together
    data_exo_pca = cat(3,data_exo_pca_day1,data_exo_pca_day2);
    data_exo_pca = cat(3,data_exo_pca,data_exo_pca_day3);
    exo_std = nanstd(data_exo_pca ,[],3);
end
%%
for tI = 1:length(tasks)% tasks = {'cd','nc','itc','lt','smb','gng','rdm'};

    % Set task to visualize
    task = tasks{tI};

    cmap = [141,211,199; % cd
        255,192,0; %255,255,179; % nc
        190,186,218; % itc
        251,128,114; % lt
        128,177,211; % smb
        253,180,98; % rdm
        179,222,105]./255; % gng
    cmap = array2table(cmap','VariableNames',{'cd','nc','itc','lt','smb','rdm','gng'}); % gng

    % NOTE: The order of the parameters here must match the order within the Stan code.
    params.itc = {'k','beta'};
    params.nc = {'weber'};
    params.lt = {'rho','beta'};
    params.smb = {'wV','wsignV_over_TU','wRU'};
    params.cd = {'criterion','sigma'};
    params.gng= {'xi','ep','b','pi','rho_rew_pun','rho_neut'};
    params.rdm = {'alpha','delta'};%{'alpha','delta','tau'};


    %%
    %     sd_frac_lrn_all = [];
    for pI = 1:length(params.(task))
        file_name = [task,'_',params.(task){pI},'.mat'];
        tmp = load(fullfile(in_dir,file_name));
        sd_frac_lrn_all.(task).(params.(task){pI}) = tmp.sd_frac_lrn;
%         x = [tmp.sd_frac_noise,tmp.sd_frac_lrn,tmp.sd_frac_exo];
        x = [tmp.sd_frac_noise,tmp.sd_frac_lrn,tmp.sd_frac_exo1,tmp.sd_frac_exo2];
%         x = [tmp.sd_frac_noise,tmp.sd_frac_lrn,tmp.sd_frac_exo1,tmp.sd_frac_exo2,tmp.sd_frac_exo3];
        figure('color','w')
% idx = isoutlier(x5,'median'); % Remove outliers based on 3MAD from the median, instead of the default 1.5IQR implemented in boxchart.
        h = boxplot(x,'jitter',0.8,'outliersize',6);
        h = boxchart(x,'Notch','off','BoxFaceColor',cmap.(task),'MarkerStyle','.','MarkerSize',20,'MarkerColor','k','BoxFaceAlpha',0.5,'BoxEdgeColor','k','BoxWidth',0.8,'LineWidth',2,'JitterOutliers','on');
        hold all
%          swarmchart(repmat(1:5,[size(x,1),1]),x,24,[1,1,1]*0.2,"filled",'MarkerFaceAlpha',0.7,'XJitterWidth',0.5,'XJitter','rand')
%          swarmchart(repmat(1:5,[size(x,1),1]),x,24,[1,1,1]*0.2,"filled",'MarkerFaceAlpha',0.7)
%          swarmchart(repmat(1:4,[size(x,1),1]),x,24,[1,1,1]*0.2)
        ylim([0,1])
        iqr_vals = prctile(x,[25,75]);
        if abs(tmp.pd_pow_a_frac)==1
            text(2,1,num2str(-tmp.pd_pow_a_frac*100),'FontSize',14) % We multiply by -1 becase a positive value means ecreaing effect on the parameter. See the powerlaw formulation to see why.
        elseif abs(tmp.pd_pow_a_frac)>=0.95 || write_all_pd
            text(2,1,num2str(-tmp.pd_pow_a_frac*100,2),'FontSize',14)
        end
        if abs(tmp.pd_A_frac_1)>=0.95 || write_all_pd
            text(3,1,['',num2str(tmp.pd_A_frac_1*100,2)],'FontSize',14)
        end
    if abs(tmp.pd_A_frac_2)>=0.95 || write_all_pd
            text(4,0.975,['', num2str(tmp.pd_A_frac_2*100,2)],'FontSize',14)
        end
%         if abs(tmp.pd_A_frac_3)>=0.95
%             text(5,0.95,['',num2str(tmp.pd_A_frac_3*100,2)],'FontSize',14)
%         end
            if add_iqr
            for ii = 1:4
                text(ii-0.1,0.95,[num2str(median(x(:,ii)),2)],'FontSize',14) % We multiply by -1 becase a positive value means ecreaing effect on the parameter. See the powerlaw formulation to see why.
                text(ii-0.2,0.93,[num2str(iqr_vals(1,ii),2),'-',num2str(iqr_vals(2,ii),2)],'FontSize',14) % We multiply by -1 becase a positive value means ecreaing effect on the parameter. See the powerlaw formulation to see why.
            end
            end

%         set(gca,'FontSize',24,'TickDir','out','xticklabel',{'Noise','Practice','Life Events'},'linewidth',2)
%         set(gca,'FontSize',24,'TickDir','out','xticklabel',{'Noise','Practice','Valence','Arousal','Self Control'},'linewidth',2)
        set(gca,'FontSize',24,'TickDir','out','xticklabel',{'Noise','Practice','Valence','Arousal'},'linewidth',2)
        title(params.(task){pI})
        axis square
        if export_figures
            set(gcf,'position',get(0,'ScreenSize'));
            fig_file = fullfile(fig_dir,sprintf('%s_relative_STD_plot_%s.png',task,params.(task){pI}));
            export_fig(gcf,fig_file,'-dpng','-r200');
            fig_file = fullfile(fig_dir,'pdf',sprintf('%s_relative_STD_plot_%s.pdf',task,params.(task){pI}));
            exportgraphics(gcf,fig_file,'ContentType','vector')
        end
%         close(gcf)
        if examine_exo_variability
            % Calculate correlation between learning contribution and variability of exo
            exo_names = {'Valence','Arousal','Self control'};
            exo_names = {'Valence','Arousal'};
            figure('color','w')
            for eI = 1:3
                subplot(1,3,eI)
                y = exo_std(eI,:);
                if strcmp(task,'lt')
                    idx = [1,4,7,8,9,11,14,15,19,23,24,25,28,32,35,36,37,47,50,53,58,60,61,65,66,68,70,74,77,83,84,88];
                    y(idx) = [];
                elseif strcmp(task,'gng')
                    idx = [1, 2, 3, 6, 8, 9, 15, 16, 17, 28, 33, 34, 37, 39, 44, 59, 60, 68, 70, 71, 74, 76, 78, 86];
                    y(idx) = [];
                end
                x = mean(sd_frac_lrn_all.(task).(params.(task){pI}),2);
                scatter(x,y,128,'filled','MarkerFaceAlpha',0.3)
                axis square
                ylabel([exo_names{eI}, ' STD'])
                xlabel('Fractional SD of Learning Terms')
                set(gca,'fontsize',14,'tickdir','out')
                title({[task,' ',params.(task){pI}];['r=',num2str(r)]})
                r = corr(x,y','Type','Pearson');
            end
        end
    end
    
continue


    
%% TEST
% % % if strcmp(task,'gng')
% % % %     continue
% % % end
% % % y = [];
% % % x = [];
% % % % sepin = 0.25;%0.4;
% % % % sepout = 1.2;
% % % % boxwidth = 0.21;%0.3;
% % % sepin = 2;%0.4;
% % % sepout = sepin*4;
% % % boxwidth = sepin*0.75;%0.3;
% % % for pI = 1:length(params.(task))
% % %     file_name = [task,'_',params.(task){pI},'.mat'];
% % %     tmp = load(fullfile(in_dir,file_name));
% % %     y = vertcat(y,[tmp.sd_frac_noise; tmp.sd_frac_lrn; tmp.sd_frac_exo]);
% % %     x = vertcat(x,x_prev+[ones(size(tmp.sd_frac_noise))-sepin; ones(size(tmp.sd_frac_noise)); ones(size(tmp.sd_frac_noise))+sepin]+(pI-1)*sepout);
% % % 
% % % params_all{end+1} = params.(task){pI};
% % % end
% % % figure(245)
% % % hold all
% % % h = boxchart(x,y,'Notch','off','BoxFaceColor',cmap.(task),'MarkerStyle','.','MarkerSize',20,'MarkerColor','k','BoxFaceAlpha',0.5,'BoxEdgeColor','k','BoxWidth',boxwidth,'LineWidth',2);
% % % xunique = unique(h.XData);
% % % 
% % % xticks_all = [xticks_all;xunique(2:3:end)];
% % % x_prev = x(end) + sepout;
% % % set(gca,'FontSize',24,'TickDir','out','xticklabel',{'Noise','Practice','Life Events'},'linewidth',2)
% % % % xunique = unique(h.XData);
% % % % set(gca,'xtick',xunique(2:3:end),'xticklabel',params.(task));
% % % title(task)
% % % % axis square
% % % % % % continue
% % % %%
% % % set(gcf,'Color','w')
% % % xlims=xlim; xlim([xlims(1)-sepin,xlims(end)+sepin]); ylim([-0.01,1])
% % % set(gca,'xtick',xticks_all,'xticklabel',params_all);
% % % 
% % %   if export_figures
% % %             set(gcf,'position',get(0,'ScreenSize'));
% % %              fig_file = fullfile(fig_dir,sprintf('%s_relative_STD_plot.png','all'));
% % %              export_fig(gcf,fig_file,'-dpng','-r200');
% % %              fig_file = fullfile(fig_dir,'pdf',sprintf('%s_relative_STD_plot.pdf','all'));
% % %              exportgraphics(gcf,fig_file,'ContentType','vector')
% % %         end
% % % % https://www.mathworks.com/matlabcentral/answers/311820-save-a-figure-as-pdf
% % % % better yet:
% % % % https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
% % % 
% % % 
% % % 
% % % % continue
% % % %% Repeat, but group each parameter separately
% % % y = [];
% % % x = [];
% % % sepin = 0.25;%0.4;
% % % sepout = 1.2;
% % % boxwidth = 0.21;%0.3;
% % % for pI = 1:length(params.(task))
% % %     file_name = [task,'_',params.(task){pI},'.mat'];
% % %     tmp = load(fullfile(in_dir,file_name));
% % %     y = vertcat(y,[tmp.sd_frac_noise; tmp.sd_frac_lrn; tmp.sd_frac_exo]);
% % %     x = vertcat(x,[ones(size(tmp.sd_frac_noise))-sepin; ones(size(tmp.sd_frac_noise)); ones(size(tmp.sd_frac_noise))+sepin]+(pI-1)*sepout);
% % % end
% % % figure('color','w')
% % % h = boxchart(x,y,'Notch','off','BoxFaceColor',cmap.(task),'MarkerStyle','.','MarkerSize',20,'MarkerColor','k','BoxFaceAlpha',0.5,'BoxEdgeColor','k','BoxWidth',boxwidth,'LineWidth',2);
% % % set(gca,'FontSize',24,'TickDir','out','xticklabel',{'Noise','Practice','Life Events'},'linewidth',2)
% % % xunique = unique(h.XData);
% % % set(gca,'xtick',xunique(2:3:end),'xticklabel',params.(task));
% % % title(task)
% % % axis square
% % % 
% % %         if export_figures_whole_task
% % %             set(gcf,'position',get(0,'ScreenSize'));
% % % %             fig_file = fullfile(fig_dir,sprintf('%s_relative_STD_plot.png',task));
% % % %             export_fig(gcf,fig_file,'-dpng','-r200');
% % % %             fig_file = fullfile(fig_dir,'pdf',sprintf('%s_relative_STD_plot.pdf',task));
% % % %             exportgraphics(gcf,fig_file,'ContentType','vector')
% % % 
% % %             set(gca,'xlim',[0,4.5])
% % %             fig_file = fullfile(fig_dir,sprintf('%s_relative_STD_plot_xlim.png',task));
% % %             export_fig(gcf,fig_file,'-dpng','-r200');
% % %             fig_file = fullfile(fig_dir,'pdf',sprintf('%s_relative_STD_plot_xlim.pdf',task));
% % %             exportgraphics(gcf,fig_file,'ContentType','vector')
% % % 
% % %         end
% % % 
% % % % % % continue

    %% Compare with accuracy
    switch task
        case 'gng'
            [accuracy,accuracy_per_week_and_condition, accuracy_per_week_and_block,accuracy_per_week_and_block_second_half,rt_per_week_and_condition] = comppheno_get_gng_accuracy(); % ~3 minutes

            exclude_subjects_with_accuracy_below_55 = true;
            if exclude_subjects_with_accuracy_below_55
                idx = find(nanmean(accuracy')<0.55);
                accuracy(idx,:) = [];
                accuracy_per_week_and_condition(idx,:,:) = [];
                accuracy_per_week_and_block(idx,:,:,:) = [];
                accuracy_per_week_and_block_second_half(idx,:,:,:) = [];
                rt_per_week_and_condition(idx,:,:) = [];
            end
        case 'rdm'
            [accuracy_per_week_and_condition,rt_per_week_and_condition] = comppheno_get_rdm_accuracy();
        case 'lt'
            [accuracy_per_week_and_condition,rt_per_week_and_cond] = comppheno_get_lt_accuracy();
            exclude_subject_choosing_safe_more_than_80_percent = true;
            if exclude_subject_choosing_safe_more_than_80_percent
                idx = [1,4,7,8,9,11,14,15,19,23,24,25,28,32,35,36,37,47,50,53,58,60,61,65,66,68,70,74,77,83,84,88];
                accuracy_per_week_and_condition(idx,:,:) = [];
                rt_per_week_and_cond(idx,:,:) = [];
            end

        case 'smb'
            [accuracy_per_week_and_condition,rt_per_week_and_condition,reward_per_week_and_condition,accuracy_per_week_and_condition_per_half] = comppheno_get_smb_accuracy();
        case 'cd'
            [accuracy_per_week_and_condition, rt_per_week_and_condition,rt_per_week_and_condition_prctile,rt_per_week_and_condition_max,rt_per_week_and_condition_outlier] = comppheno_get_cd_accuracy();
        case 'nc'
            [accuracy_per_week,rt_per_week] = comppheno_get_nc_accuracy();
            accuracy_per_week_and_condition = accuracy_per_week;
        case 'itc'
            [rt_per_week] = comppheno_get_itc_accuracy();
            accuracy_per_week_and_condition = rt_per_week;
        otherwise
            warning('Please enter a valid task name.')
    end


    %%
    % accuracy_per_week = nanmean(accuracy_per_week_and_condition,3);
    n_conditions = size(accuracy_per_week_and_condition,3);

    figure('color','w')
    for cI = 1:n_conditions
        accuracy_per_week = squeeze(accuracy_per_week_and_condition(:,:,cI));
        accuracy_diff = [];
        for sI = 1:size(accuracy_per_week,1)
            idx = find(accuracy_per_week(sI,:)>0);

            accuracy_diff(sI) = accuracy_per_week(sI,idx(end)) -  accuracy_per_week(sI,idx(1));
        end

        figure('color','w')
        for pI = 1:length(params.(task))
            subplot(1,length(params.(task)),pI)
            %         subplot(4,length(params.(task)),pI+5*(cI-1))

            x = mean(sd_frac_lrn_all.(task).(params.(task){pI}),2);
            y = accuracy_diff;
            scatter(x,y,128,'filled','MarkerFaceAlpha',0.3)
            axis square
            if pI==1
                ylabel('\DeltaAccuracy')
                xlabel('Fractional SD of Learning Terms')
            end
            set(gca,'fontsize',14,'tickdir','out')
            r = corr(x,y','Type','Pearson');
            pval = perm_test(x,y',r,5000);

            title({[task,' ',params.(task){pI}]; ['r_P=',num2str(r,2), ' p_p_e_r_m=',num2str(pval,2)]})
        end
        if export_learning_figures
            set(gcf,'position',get(0,'ScreenSize'));
            fig_file = fullfile(fig_dir,'relative_learning_STD_vs_accuracyDiff',sprintf('%s_relative_learning_STD_vs_accuracyDiff_%s.png',task,num2str(cI)));
            export_fig(gcf,fig_file,'-dpng','-r200');
            fig_file = fullfile(fig_dir,'relative_learning_STD_vs_accuracyDiff','pdf',sprintf('%s_relative_learning_STD_vs_accuracyDiff_%s.pdf',task,num2str(cI)));
            exportgraphics(gcf,fig_file,'ContentType','vector')
        end
    end

    % close all

end




%% Demonstrate the time series of different termsedit 

task = 'gng';
pI = 3;
[~,sI] = max(tmp.sd_frac_lrn);

file_name = [task,'_',params.(task){pI},'.mat'];
tmp = load(fullfile(in_dir,file_name));

figure('color','w')
for sI = 1:49
subplot(7,7,sI)
hold all
plot(tmp.noise_mean(sI,:))
plot(tmp.lrn_mean(sI,:))
plot(tmp.exo1_mean(sI,:))
plot(tmp.exo2_mean(sI,:))
axis square
box off
% set(gca,'fontsize',14,'tickdir','out')
title(num2str(sI))
end


figure('color','w')
[~,sI] = max(tmp.sd_frac_exo2);
hold all
plot(tmp.noise_mean(sI,:))
plot(tmp.lrn_mean(sI,:))
plot(tmp.exo1_mean(sI,:))
plot(tmp.exo2_mean(sI,:))
axis square
box off
set(gca,'fontsize',14,'tickdir','out')
title(num2str(sI))
