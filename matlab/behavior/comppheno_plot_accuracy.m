
addpath('/n/home10/rschurr/Matlab_code/external_toolboxes/altmany-export_fig-3e1e29a'); % For exporting figures

matplotlib_colormaps

plot_rt = false; % True: plot RT; False: plot accuracy

add_power_law_fit = true;
n_perms = 10000; % Number of random permutations of time to calculate a p-value for the adjusted R^2 of the power-law fit
set_ylim = false; % Set specific ylim to the figures
set_ylim_tight = true;
ylims_all = [0.5,1];

% General 3-term power law settings
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display='off';
opts.Upper = [Inf 0 inf];

%%
if plot_rt
    label = 'RT [ms]';
else
    label = 'accuracy';
end

% set(gcf,'position',get(0,'screensize'))

export_figures = true;
fig_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/comp_pheno_figs/accuracy_rt_over_time/';
sub_dir = '';
if ~ set_ylim
    sub_dir = 'variable_ylim';
end
if ~exist(fullfile(fig_dir,sub_dir),'dir')
    mkdir(fullfile(fig_dir,sub_dir))
end
if ~exist(fullfile(fig_dir,sub_dir,'pdf'),'dir')
    mkdir(fullfile(fig_dir,sub_dir,'pdf'))
end



%% Plot RDM accuracy
[rdm_accuracy_per_week_and_coherence,rdm_rt_per_week_and_coherence] = comppheno_get_rdm_accuracy();

%% RDM heat map
if plot_rt
    tmp = mean(rdm_rt_per_week_and_coherence,3); % Average across coherence values
else
    tmp = mean(rdm_accuracy_per_week_and_coherence,3); % Average across coherence values
end
% figure('color','w')
% nanimagesc(tmp)
% caxis(prctile(tmp(:),[2,98]))
% colorbar
% set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
% xlabel('session')
% ylabel('subject')
% title(['RDM ' label])
% axis square
%
% set(gcf,'position',get(0,'screensize'))
% if export_figures
%     fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_RDM_hetamp_',label,'.png']);
%     export_fig(gcf,fig_file,'-dpng','-r200');
% end


% RDM curve
figure('color','w')
% h = plot(nanmean(tmp),'.','linew',2);
h = scatter([1:12],nanmean(tmp),128,[0 0 0],'filled');
hold all
herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h.CData,'LineStyle','none','LineWidth',2);
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel(label)
title(['RDM ' label])
axis square
% h = plot(tmp','linew',0.5,'Color',0.5*[1,1,1]); % Add individual subjects
pause(1)
if add_power_law_fit
    x = 1:12;
    y = nanmean(tmp);
    [xData, yData] = prepareCurveData( x, y );


    [fitresult, gof_orig] = fit( xData, yData, ft, opts );
    t = x(1):0.1:x(end);
    hfit = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');

    hfit.Color = h.CData;%0.5*[1,1,1]
    hfit.LineWidth = 2;
    subtitle(['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)])
 
    % Permute time and repeat the power law fit, to get a p-value for the
    % adjusted R^2 value
    rsquare_perms = nan(1,n_perms);
    adjrsquare_perms = nan(1,n_perms);
    for permI = 1:n_perms
        perm = randperm(12);
        y_perm = y(perm);
        [xData, yData] = prepareCurveData( x, y_perm);
        opts.Upper = [Inf 0 inf];
        [fitresult, gof] = fit( xData, yData, ft, opts ); 
        rsquare_perms(permI) = gof.rsquare;
        adjrsquare_perms(permI) = gof.adjrsquare;
    end
    prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
    prctile_str = num2str(prctile_num,2);
    subtitle({  ['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]; ['p=', num2str((100-prctile_num)/100)]  })
    
end

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_RDM_curve_',label,'.png']);
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    if set_ylim_tight
        ylim([0.76,0.84]);
    end
    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end

pause(0.5)
clear h
clear hfit

% RDM curves (one for each coherence value)
figure('color','w')
cmap = copper(size(rdm_accuracy_per_week_and_coherence,3));
for coh = 1:size(rdm_accuracy_per_week_and_coherence,3)

    if plot_rt
        tmp = rdm_rt_per_week_and_coherence(:,:,coh);
    else
        tmp = rdm_accuracy_per_week_and_coherence(:,:,coh);
    end

    %     h(coh) = plot(nanmean(tmp),'.','linew',2,'Color',cmap(coh,:));
    h(coh) = scatter([1:12],nanmean(tmp),128,cmap(coh,:),'filled');

    hold on
    errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h(coh).CData,'LineStyle','none','LineWidth',2)
    if add_power_law_fit
        x = 1:12;
        y = nanmean(tmp);
        [xData, yData] = prepareCurveData( x, y );

        [fitresult, gof_orig] = fit( xData, yData, ft, opts );
        t = x(1):0.1:x(end);
        hfit(coh) = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
        radj(coh) = gof_orig.adjrsquare;
        hfit(coh).Color = h(coh).CData;%0.5*[1,1,1]
        hfit(coh).LineWidth = 2;

    % Permute time and repeat the power law fit, to get a p-value for the
    % adjusted R^2 value
    rsquare_perms = nan(1,n_perms);
    adjrsquare_perms = nan(1,n_perms);
    for permI = 1:n_perms
        perm = randperm(12);
        y_perm = y(perm);
        [xData, yData] = prepareCurveData( x, y_perm);
        opts.Upper = [Inf 0 inf];
        [fitresult, gof] = fit( xData, yData, ft, opts ); 
        rsquare_perms(permI) = gof.rsquare;
        adjrsquare_perms(permI) = gof.adjrsquare;
    end
    prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
    pval_str{coh}=['p=', num2str((100-prctile_num)/100)];
    
    end
end
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel(label)
title(['RDM ' label])
axis square
% legend([h(1),h(2),h(3),h(4)],{'5','10','35','50'})
% legend([hfit],{['5 R^2_a_d_j=' num2str(round(radj(1)*10)/10,2),' ',pval_str{1}],...
%     ['10 R^2_a_d_j='  num2str(round(radj(2)*10)/10,2),' ',pval_str{2}],...
%     ['35 R^2_a_d_j=' num2str(round(radj(3)*10)/10,2),' ',pval_str{3}],...
%     ['50 R^2_a_d_j=' num2str(round(radj(4)*10)/10,2),' ',pval_str{4}]})

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_RDM_curves_',label,'.png']);
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end


%%
%% Plot CD accuracy
if ~plot_rt % CD does not have RT data to plot
    [cd_accuracy_per_week_and_condition] = comppheno_get_cd_accuracy(); % No RT data for the CD task

    %% CD heat map
    tmp = mean(cd_accuracy_per_week_and_condition,3); % Average across set size and number of targets

%     figure('color','w')
%     nanimagesc(tmp)
%     caxis(prctile(tmp(:),[2,98]))
%     colorbar
%     set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
%     xlabel('session')
%     ylabel('subject')
%     title(['CD ' label])
%     axis square
% 
%     set(gcf,'position',get(0,'screensize'))
%     if export_figures
%         fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_CD_hetamp_',label,'.png']);
%         if set_ylim
%             fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
%             ylim(ylims_all);
%         end
%         if set_ylim_tight
%             ylim([0.73,0.79])
%         end
%         export_fig(gcf,fig_file,'-dpng','-r200');
%     end

    % CD curve
    figure('color','w')
    %     h = plot(nanmean(tmp),'.','linew',2);
    h = scatter([1:12],nanmean(tmp),128,[0,0,0],'filled');
    hold all
    herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h.CData,'LineStyle','none','LineWidth',2)
    ylims = ylim
    ylims(1) = floor(ylims(1)*10)/10;
    ylims(2) = ceil(ylims(2)*10)/10;
    ylim(ylims); % Round accuracy values to closest 10,20, etc.
    set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
    xlabel('session')
    ylabel(label)
    title(['CD ' label])
    axis square

    if add_power_law_fit
        x = 1:12;
        y = nanmean(tmp);
        [xData, yData] = prepareCurveData( x, y );
        % Fit model to data.
        [fitresult, gof_orig] = fit( xData, yData, ft, opts );
        %     hfit = plot(fitresult);
        t = 1:0.1:12;
        hfit = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
        hfit.Color = h.CData;%0.5*[1,1,1]
        hfit.LineWidth = 2;
        subtitle(['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)])

        % Permute time and repeat the power law fit, to get a p-value for the
        % adjusted R^2 value
        rsquare_perms = nan(1,n_perms);
        adjrsquare_perms = nan(1,n_perms);
        for permI = 1:n_perms
            perm = randperm(12);
            y_perm = y(perm);
            [xData, yData] = prepareCurveData( x, y_perm);
            opts.Upper = [Inf 0 inf];
            [fitresult, gof] = fit( xData, yData, ft, opts );
            rsquare_perms(permI) = gof.rsquare;
            adjrsquare_perms(permI) = gof.adjrsquare;
        end
        prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
        prctile_str = num2str(prctile_num,2);
        subtitle({  ['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]; ['p=', num2str((100-prctile_num)/100)]  })
    end
figure(4)
    set(gcf,'position',get(0,'screensize'))
    if export_figures
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_CD_curve_',label,'.png']);
        if set_ylim
            fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
            ylim(ylims_all);
        end
        if set_ylim_tight
             ylim([0.73,0.79])
         end

        export_fig(gcf,fig_file,'-dpng','-r200');
        fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

    end

    pause(0.5)
    clear h
    clear hfit

    % CD curves (one for each condition)
    figure('color','w')
    cmap = copper(size(cd_accuracy_per_week_and_condition,3));
    cmap(end,:) = cmap(end-1,:); % Make both 8-items the same color
    for cond = 1:size(cd_accuracy_per_week_and_condition,3)
        if plot_rt
            tmp = cd_rt_per_week_and_condition(:,:,cond);
        else
            tmp = cd_accuracy_per_week_and_condition(:,:,cond);
        end
        %h(cond) = plot(nanmean(tmp),'.','linew',2,'Color',cmap(cond,:));
        h(cond) = scatter([1:12],nanmean(tmp),128,cmap(cond,:),'filled');
        hold on
        errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'.','Color',h(cond).CData)

        if add_power_law_fit
            x = 1:12;
            y = nanmean(tmp);
            [xData, yData] = prepareCurveData( x, y );
            % Fit model to data.
            [fitresult, gof_orig] = fit( xData, yData, ft, opts );
            %     hfit = plot(fitresult);
            radj(cond) = gof_orig.adjrsquare;
            t = 1:0.1:12;
            hfit(cond) = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
            hfit(cond).Color = h(cond).CData;%0.5*[1,1,1]
            hfit(cond).LineWidth = 2;

            % Permute time and repeat the power law fit, to get a p-value for the
            % adjusted R^2 value
            rsquare_perms = nan(1,n_perms);
            adjrsquare_perms = nan(1,n_perms);
            for permI = 1:n_perms
                perm = randperm(12);
                y_perm = y(perm);
                [xData, yData] = prepareCurveData( x, y_perm);
                opts.Upper = [Inf 0 inf];
                [fitresult, gof] = fit( xData, yData, ft, opts ); 
                rsquare_perms(permI) = gof.rsquare;
                adjrsquare_perms(permI) = gof.adjrsquare;
            end
            prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
            pval_str{cond}=['p=', num2str((100-prctile_num)/100)];
        end

    end
    ylims = ylim;
    ylims(1) = floor(ylims(1)*10)/10;
    ylims(2) = ceil(ylims(2)*10)/10;
    ylim(ylims); % Round accuracy values to closest 10,20, etc.
    set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
    xlabel('session')
    ylabel(label)
    title(['CD ' label])
    axis square
%     legend([h(1),h(2),h(3),h(4),h(5)],{'3','4','6','8','8vary'})
    % h(5).LineStyle = '--';
    hfit(5).LineStyle = '--';

%     legend([hfit],{['3 R^2_a_d_j=' num2str(round(radj(1)*10)/10,2),' ',pval_str{1}],...
%         ['4 R^2_a_d_j='  num2str(round(radj(2)*10)/10,2),' ',pval_str{2}],...
%         ['6 R^2_a_d_j=' num2str(round(radj(3)*10)/10,2),' ',pval_str{3}],...
%         ['8 R^2_a_d_j=' num2str(round(radj(4)*10)/10,2),' ',pval_str{4}],...
%         ['8_v_a_r_y R^2_a_d_j=' num2str(round(radj(5)*10)/10,2),' ',pval_str{5}]})

    set(gcf,'position',get(0,'screensize'))
    if export_figures
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_CD_curves_',label,'.png']);
        if set_ylim
            fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
            ylim(ylims_all);
        end

        export_fig(gcf,fig_file,'-dpng','-r200');
        
    end

end
%%
%% Plot NC accuracy
[nc_accuracy_per_week,nc_rt_per_week] = comppheno_get_nc_accuracy();

%% NC heat map

if plot_rt
    tmp = nc_rt_per_week;
else
    tmp = nc_accuracy_per_week;
end

% figure('color','w')
% nanimagesc(tmp)
% caxis(prctile(tmp(:),[2,98]))
% colorbar
% set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
% xlabel('session')
% ylabel('subject')
% title('NC Accuracy')
% axis square
% 
% set(gcf,'position',get(0,'screensize'))
% if export_figures
%     fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_NC_hetamp_',label,'.png']);
%     if set_ylim
%         fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
%         ylim(ylims_all);
%     end
%     export_fig(gcf,fig_file,'-dpng','-r200');
%     fig_file = strrep(fig_file,'.png','.pdf');
% %     fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
%     exportgraphics(gcf,fig_file,'ContentType','vector')
% 
% end


clear h
clear hfit

% NC curve
figure('color','w')
% h = plot(nanmean(tmp),'.','linew',2);
h  = scatter([1:12],nanmean(tmp),128,[0,0,0],'filled');

hold all
herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h.CData,'LineStyle','none','LineWidth',2);
ylims = ylim
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel(label)
title('NC accuracy')
axis square
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel(label)
title(['NC ' label])
axis square


if add_power_law_fit
    x = 1:12;
    y = nanmean(tmp);
    [xData, yData] = prepareCurveData( x, y );
    % Fit model to data.
    [fitresult, gof_orig] = fit( xData, yData, ft, opts );
    %     hfit = plot(fitresult);
    t = 1:0.1:12;
    hfit = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
    hfit.Color = h.CData;%0.5*[1,1,1]
    hfit.LineWidth = 2;
    subtitle(['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]);

    % Permute time and repeat the power law fit, to get a p-value for the
    % adjusted R^2 value
    rsquare_perms = nan(1,n_perms);
    adjrsquare_perms = nan(1,n_perms);
    for permI = 1:n_perms
        perm = randperm(12);
        y_perm = y(perm);
        [xData, yData] = prepareCurveData( x, y_perm);
        opts.Upper = [Inf 0 inf];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare_perms(permI) = gof.rsquare;
        adjrsquare_perms(permI) = gof.adjrsquare;
    end
    prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
    prctile_str = num2str(prctile_num,2);
    subtitle({  ['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]; ['p=', num2str((100-prctile_num)/100)]  })
end

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_NC_curve_',label,'.png']);
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
         if set_ylim_tight
         ylim([0.83,0.87])
     end

    export_fig(gcf,fig_file,'-dpng','-r200');
        fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end

%%
%% Plot GNG accuracy
[gng_accuracy,gng_accuracy_per_week_and_condition, gng_accuracy_per_week_and_block,gng_accuracy_per_week_and_block_second_half,gng_rt_per_week_and_condition] = comppheno_get_gng_accuracy(); % ~3 minutes

exclude_subjects_with_accuracy_below_55 = true;
if exclude_subjects_with_accuracy_below_55
    idx = find(nanmean(gng_accuracy')<0.55);
    gng_accuracy(idx,:) = [];
    gng_accuracy_per_week_and_condition(idx,:,:) = [];
    gng_accuracy_per_week_and_block(idx,:,:,:) = [];
    gng_accuracy_per_week_and_block_second_half(idx,:,:,:) = [];
    gng_rt_per_week_and_condition(idx,:,:) = [];
end


%% GNG heat map
if plot_rt
    tmp = mean(gng_rt_per_week_and_condition,3); % Average across set size and number of targets
else
    tmp = mean(gng_accuracy_per_week_and_condition,3); % Average across set size and number of targets
end
% figure('color','w')
% nanimagesc(tmp)
% caxis(prctile(tmp(:),[2,98]))
% colorbar
% set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
% xlabel('session')
% ylabel('subject')
% title(['GNG ' label])
% axis square
% 
% set(gcf,'position',get(0,'screensize'))
% if export_figures
%     fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_GNG_hetamp_',label,'.png']);
%     if exclude_subjects_with_accuracy_below_55
%         fig_file = strrep(fig_file,'.png','_exclude_subjects_with_accuracy_below_55.png');
%     end
%     if set_ylim
%         fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
%         ylim(ylims_all);
%     end
% 
% export_fig(gcf,fig_file,'-dpng','-r200');
% fig_file = strrep(fig_file,'.png','.pdf');
% % fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
% exportgraphics(gcf,fig_file,'ContentType','vector')
% end


% GNG curve
figure('color','w')
% h = plot(nanmean(tmp),'.','linew',2);
h = scatter([1:12],nanmean(tmp),128,[0,0,0],'filled');
hold all
herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h.CData,'LineStyle','none','LineWidth',2);
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel(label)
title('GNG accuracy')
axis square
if add_power_law_fit
    x = 1:12;
    y = nanmean(tmp);
    [xData, yData] = prepareCurveData( x, y );
    % Fit model to data.
    [fitresult, gof_orig] = fit( xData, yData, ft, opts );
    %     hfit = plot(fitresult);
    t = 1:0.1:12;
    hfit = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
    hfit.Color = h.CData;%0.5*[1,1,1]
    hfit.LineWidth = 2;
    subtitle(['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)])

    % Permute time and repeat the power law fit, to get a p-value for the
    % adjusted R^2 value
    rsquare_perms = nan(1,n_perms);
    adjrsquare_perms = nan(1,n_perms);
    for permI = 1:n_perms
        perm = randperm(12);
        y_perm = y(perm);
        [xData, yData] = prepareCurveData( x, y_perm);
        opts.Upper = [Inf 0 inf];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare_perms(permI) = gof.rsquare;
        adjrsquare_perms(permI) = gof.adjrsquare;
    end
    prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
    prctile_str = num2str(prctile_num,2);
    subtitle({  ['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]; ['p=', num2str((100-prctile_num)/100)]  })
end

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_GNG_curve_',label,'.png']);
    if exclude_subjects_with_accuracy_below_55
        fig_file = strrep(fig_file,'.png','_exclude_subjects_with_accuracy_below_55.png');
    end
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    export_fig(gcf,fig_file,'-dpng','-r200');

 fig_file = strrep(fig_file,'.png','.pdf');
%  fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
     exportgraphics(gcf,fig_file,'ContentType','vector')

end

pause(0.5)
clear h
clear hfit


% GNG curves (one for each condition)
figure('color','w')
cmap = copper(size(gng_accuracy_per_week_and_condition,3));
for cond = 1:size(gng_accuracy_per_week_and_condition,3)

    if plot_rt
        tmp = gng_rt_per_week_and_condition(:,:,cond);
    else
        tmp = gng_accuracy_per_week_and_condition(:,:,cond);
    end

    %     h(cond) = plot(nanmean(tmp),'.','linew',2,'Color',cmap(cond,:));
    h(cond)  = scatter([1:12],nanmean(tmp),128,cmap(cond,:),'filled');
    hold on
    errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h(cond).CData,'LineStyle','none','LineWidth',2); % Note: For RT, the NoGo conditions reflect the erroneous Go responses (because NoGo have no RT data)

    if add_power_law_fit
        x = 1:12;
        y = nanmean(tmp);

        [xData, yData] = prepareCurveData( x, y );
        [fitresult, gof_orig] = fit( xData, yData, ft, opts );
        radj(cond) = gof_orig.adjrsquare;
        t = x(1):0.1:x(end);

        hfit(cond) = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
        % % %     hfit(cond) = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
        hfit(cond) = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');

        hfit(cond).Color = h(cond).CData;%0.5*[1,1,1]
        hfit(cond).LineWidth = 2;

            % Permute time and repeat the power law fit, to get a p-value for the
            % adjusted R^2 value
            rsquare_perms = nan(1,n_perms);
            adjrsquare_perms = nan(1,n_perms);
            for permI = 1:n_perms
                perm = randperm(12);
                y_perm = y(perm);
                [xData, yData] = prepareCurveData( x, y_perm);
                opts.Upper = [Inf 0 inf];
                [fitresult, gof] = fit( xData, yData, ft, opts ); 
                rsquare_perms(permI) = gof.rsquare;
                adjrsquare_perms(permI) = gof.adjrsquare;
            end
            prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
            pval_str{cond}=['p=', num2str((100-prctile_num)/100)];
    end

end
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel(label)
title('GNG accuracy')
axis square
% legend([h(1),h(2),h(3),h(4)],{'Go2win','NoGo2win','Go2avoid','NoGo2avoid'})

% subtitle(['R^2_a_d_j: Go2win=', num2str(round(radj(1)*10)/10,2),...
%     '; NoGo2win=', num2str(round(radj(2)*10)/10,2)])

% legend(hfit,{['Go2win R^2_a_d_j=' num2str(round(radj(1)*10)/10,2),' ',pval_str{1}],...
%     ['NoGo2win R^2_a_d_j='  num2str(round(radj(2)*10)/10,2),' ',pval_str{2}],...
%     ['Go2avoid R^2_a_d_j=' num2str(round(radj(3)*10)/10,2),' ',pval_str{3}],...
%     ['NoGo2avoid R^2_a_d_j=' num2str(round(radj(4)*10)/10,2),' ',pval_str{4}]})

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_GNG_curves_',label,'.png']);
    if exclude_subjects_with_accuracy_below_55
        fig_file = strrep(fig_file,'.png','_exclude_subjects_with_accuracy_below_55.png');
    end
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end

    export_fig(gcf,fig_file,'-dpng','-r200');

 fig_file = strrep(fig_file,'.png','.pdf');
%  fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
     exportgraphics(gcf,fig_file,'ContentType','vector')

end

%% Plot SMB accuracy
[smb_accuracy_per_week_and_condition,smb_rt_per_week_and_condition,smb_reward_per_week_and_condition,smb_accuracy_per_week_and_condition_per_half] = comppheno_get_smb_accuracy();

% SMB heat map
if plot_rt
    tmp = nanmean(smb_rt_per_week_and_condition,3); % Average across set size and number of targets
else
    tmp = nanmean(smb_accuracy_per_week_and_condition,3); % Average across set size and number of targets
end

figure('color','w')
imagesc(tmp)
caxis(prctile(tmp(:),[2,98]))
colorbar
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel('subject')
title('Accuracy')
axis square

set(gcf,'position',get(0,'screensize'))
if export_figures
    if plot_rt
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_SMB_hetamp_RT [ms].png']);
    else
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_SMB_hetamp_accuracy.png']);
    end
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end

    export_fig(gcf,fig_file,'-dpng','-r200');
end

% SMB curve
figure('color','w')
% h = plot(nanmean(tmp),'.','linew',2);
h = scatter([1:12],nanmean(tmp),128,[0,0,0],'filled');
hold all
herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h.CData,'LineStyle','none','LineWidth',2);
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
if plot_rt
    ylabel('RT [ms]')
    title('Bandit RT')
else
    ylabel('accuracy')
    title('Bandit accuracy')
end
axis square

if add_power_law_fit
    x = 1:12;
    y = nanmean(tmp);
    [xData, yData] = prepareCurveData( x, y );
    % Fit model to data.
    [fitresult, gof_orig] = fit( xData, yData, ft, opts );
    %     hfit = plot(fitresult);
    t = 1:0.1:12;
    hfit = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
    hfit.Color = h.CData;%0.5*[1,1,1]
    hfit.LineWidth = 2;
    subtitle(['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)])

    % Permute time and repeat the power law fit, to get a p-value for the
    % adjusted R^2 value
    rsquare_perms = nan(1,n_perms);
    adjrsquare_perms = nan(1,n_perms);
    for permI = 1:n_perms
        perm = randperm(12);
        y_perm = y(perm);
        [xData, yData] = prepareCurveData( x, y_perm);
        opts.Upper = [Inf 0 inf];
        [fitresult, gof] = fit( xData, yData, ft, opts ); 
        rsquare_perms(permI) = gof.rsquare;
        adjrsquare_perms(permI) = gof.adjrsquare;
    end
    prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
    prctile_str = num2str(prctile_num,2);
    subtitle({  ['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]; ['p=', num2str((100-prctile_num)/100)]  })

end

set(gcf,'position',get(0,'screensize'))
if export_figures
    if plot_rt
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_SMB_curve_RT [ms].png']);
    else
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_SMB_curve_accuracy.png']);
    end
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    if set_tight_ylim
        ylim([0.76,0.82]);
    end

    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end


pause(0.5)
% SMB curves (one for each condition)
figure('color','w')
cmap = copper(size(smb_accuracy_per_week_and_condition,3));
clear h
clear hfit
for cond = 1:size(smb_accuracy_per_week_and_condition,3)
    if cond==2
        continue % Skip RS, since it combined with SR in cond=1
    end
    if plot_rt
        tmp = smb_rt_per_week_and_condition(:,:,cond);
    else
        tmp = smb_accuracy_per_week_and_condition(:,:,cond);
        if cond==1
        tmp = mean(smb_accuracy_per_week_and_condition(:,:,1:2),3);
        end
    end

    %     h(cond) = plot(nanmean(tmp),'.','linew',2,'Color',cmap(cond,:));
    h(cond) = scatter([1:12],nanmean(tmp),128,cmap(cond,:),'filled');
    hold on
    errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h(cond).CData,'LineStyle','none','LineWidth',2)

    if add_power_law_fit
        x = 1:12;
        y = nanmean(tmp);
        [xData, yData] = prepareCurveData( x, y );
        % Fit model to data.
        [fitresult, gof_orig] = fit( xData, yData, ft, opts );
        %     hfit = plot(fitresult);
        radj(cond) = gof_orig.adjrsquare;
        t = 1:0.1:12;
        hfit(cond) = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
        hfit(cond).Color = h(cond).CData;%0.5*[1,1,1]
        hfit(cond).LineWidth = 2;

            % Permute time and repeat the power law fit, to get a p-value for the
            % adjusted R^2 value
            rsquare_perms = nan(1,n_perms);
            adjrsquare_perms = nan(1,n_perms);
            for permI = 1:n_perms
                perm = randperm(12);
                y_perm = y(perm);
                [xData, yData] = prepareCurveData( x, y_perm);
                opts.Upper = [Inf 0 inf];
                [fitresult, gof] = fit( xData, yData, ft, opts ); 
                rsquare_perms(permI) = gof.rsquare;
                adjrsquare_perms(permI) = gof.adjrsquare;
            end
            prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
            pval_str{cond}=['p=', num2str((100-prctile_num)/100)];
    end

end
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel('accuracy')
title('Accuracy')
axis square
% legend([h(1),h(2),h(3),h(4)],{'RS','SR','RR','SS'})
% legend([hfit],{['RS R^2_a_d_j=' num2str(round(radj(1)*10)/10,2),' ',pval_str{1}],...
%     ['SR R^2_a_d_j='  num2str(round(radj(2)*10)/10,2),' ',pval_str{2}],...
%     ['RR R^2_a_d_j=' num2str(round(radj(3)*10)/10,2),' ',pval_str{3}],...
%     ['SS R^2_a_d_j=' num2str(round(radj(4)*10)/10,2),' ',pval_str{4}]})

set(gcf,'position',get(0,'screensize'))
if export_figures
    if plot_rt
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_SMB_curves_RT [ms].png']);
    else
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_SMB_curves_accuracy.png']);
    end
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end

    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end


%% ITC
% NOTE: For ITC there is no definition of accuracy, so we can only plot RT.[smb_accuracy_per_week_and_condition,smb_rt_per_week_and_condition,smb_reward_per_week_and_condition] = comppheno_get_smb_accuracy();

[itc_rt_per_week] = comppheno_get_itc_accuracy();

% ITC heat map
tmp = itc_rt_per_week;

% figure('color','w')
% nanimagesc(tmp)
% caxis(prctile(tmp(:),[2,98]))
% colorbar
% set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
% xlabel('session')
% ylabel('subject')
% title('ITC accuracy')
% axis square
% 
% set(gcf,'position',get(0,'screensize'))
% if export_figures
%     fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_ITC_hetamp_RT [ms].png']);
%     if set_ylim
%         fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
%         ylim(ylims_all);
%     end
%     export_fig(gcf,fig_file,'-dpng','-r200');
% end

% ITC curve
figure('color','w')
h = plot(nanmean(tmp),'.','linew',2);
h = scatter([1:12],nanmean(tmp),128,[0,0,0],'filled');
hold all
herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h.CData,'LineStyle','none','LineWidth',2);
ylims = ylim
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
ylabel('RT [ms]')
title('ITC RT')
axis square

if add_power_law_fit
    x = 1:12;
    y = nanmean(tmp);
    [xData, yData] = prepareCurveData( x, y );
    % Fit model to data.
    [fitresult, gof_orig] = fit( xData, yData, ft, opts );
    %     hfit = plot(fitresult);
    t = 1:0.1:12;
    hfit = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
    hfit.Color = h.CData;%0.5*[1,1,1]
    hfit.LineWidth = 2;
    subtitle(['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)])

    % Permute time and repeat the power law fit, to get a p-value for the
    % adjusted R^2 value
    rsquare_perms = nan(1,n_perms);
    adjrsquare_perms = nan(1,n_perms);
    for permI = 1:n_perms
        perm = randperm(12);
        y_perm = y(perm);
        [xData, yData] = prepareCurveData( x, y_perm);
        opts.Upper = [Inf 0 inf];
        [fitresult, gof] = fit( xData, yData, ft, opts ); 
        rsquare_perms(permI) = gof.rsquare;
        adjrsquare_perms(permI) = gof.adjrsquare;
    end
    prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
    prctile_str = num2str(prctile_num,2);
    subtitle({  ['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]; ['p=', num2str((100-prctile_num)/100)]  })
end

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_ITC_curve_RT [ms].png']);
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end

%% Plot LT "accuracy"
[lt_accuracy_per_week_and_cond,lt_rt_per_week_and_cond] = comppheno_get_lt_accuracy();

exclude_subject_choosing_safe_more_than_80_percent = true;
if exclude_subject_choosing_safe_more_than_80_percent 
    idx = [1,4,7,8,9,11,14,15,19,23,24,25,28,32,35,36,37,47,50,53,58,60,61,65,66,68,70,74,77,83,84,88];
    lt_accuracy_per_week_and_cond(idx,:,:) = [];
    lt_rt_per_week_and_cond(idx,:,:) = [];
end

% LT heat map
if plot_rt
    tmp = nanmean(lt_rt_per_week_and_cond,3); % Average across set size and number of targets
else
    tmp = nanmean(lt_accuracy_per_week_and_cond,3); % Average across size condition

end

% figure('color','w')
% nanimagesc(tmp)
% caxis(prctile(tmp(:),[2,98]))
% colorbar
% set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
% xlabel('session')
% ylabel('subject')
% title('Accuracy')
% axis square
% 
% set(gcf,'position',get(0,'screensize'))
% if export_figures
%     if plot_rt
%         fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_LT_hetamp_RT [ms].png']);
%     else
%         fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_LT_hetamp_accuracy.png']);
%     end
% 
% 
%     export_fig(gcf,fig_file,'-dpng','-r200');
% end

% LT curve
figure('color','w')
% h = plot(nanmean(tmp),'.k','linew',2);
h = scatter(1:12,nanmean(tmp),128,[0 0 0],'filled');
hold all
herr = errorbar(1:12,nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h.CData,'LineStyle','none','LineWidth',2);
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
if plot_rt
    ylabel('RT [ms]')
    title('LT RT')
else
    ylabel('accuracy')
    title('LT accuracy')
end
axis square

if add_power_law_fit
    x = 1:12;
    y = nanmean(tmp);
    [xData, yData] = prepareCurveData( x, y );
    % Fit model to data.
    [fitresult, gof_orig] = fit( xData, yData, ft, opts );
    %     hfit = plot(fitresult);
    t = 1:0.1:12;
    hfit = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
    hfit.Color = h.CData;%0.5*[1,1,1]
    hfit.LineWidth = 2;
    subtitle(['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)])

    % Permute time and repeat the power law fit, to get a p-value for the
    % adjusted R^2 value
    rsquare_perms = nan(1,n_perms);
    adjrsquare_perms = nan(1,n_perms);
    for permI = 1:n_perms
        perm = randperm(12);
        y_perm = y(perm);
        [xData, yData] = prepareCurveData( x, y_perm);
        opts.Upper = [Inf 0 inf];
        [fitresult, gof] = fit( xData, yData, ft, opts ); 
        rsquare_perms(permI) = gof.rsquare;
        adjrsquare_perms(permI) = gof.adjrsquare;
    end
    prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
    prctile_str = num2str(prctile_num,2);
    subtitle({  ['R^2_a_d_j = ', num2str(gof_orig.adjrsquare,2)]; ['p=', num2str((100-prctile_num)/100)]  })
end

set(gcf,'position',get(0,'screensize'))
if export_figures
    if plot_rt
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_LT_curve_RT [ms].png']);
    else
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_LT_curve_accuracy.png']);
    end
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    if set_ylim_tight
        ylim([0.7 0.82]);
    end

    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end


pause(0.5)
%% LT curves (one for each condition)
jittervals = [-0.1,0,0.1];
figure('color','w')
cmap = copper(size(lt_accuracy_per_week_and_cond,3));
clear h
clear hfit
for cond = 1:size(lt_accuracy_per_week_and_cond,3)
    if plot_rt
        tmp = lt_rt_per_week_and_cond(:,:,cond);
    else
        tmp = lt_accuracy_per_week_and_cond(:,:,cond);
    end

    %     h(cond) = plot([1:12]+jittervals(cond),nanmean(tmp),'.','linew',2,'Color',cmap(cond,:));
    h(cond) = scatter([1:12]+jittervals(cond),nanmean(tmp),128,cmap(cond,:),'filled');

    hold on
    errorbar([1:12]+jittervals(cond),nanmean(tmp),nanstd(tmp)./sqrt(90),'Color',h(cond).CData,'LineStyle','none','LineWidth',2)

    if add_power_law_fit
        x = 1:12;
        y = nanmean(tmp);
        [xData, yData] = prepareCurveData( x, y );
        % Fit model to data.
        [fitresult, gof_orig] = fit( xData, yData, ft, opts );
        %     hfit = plot(fitresult);
        radj(cond) = gof_orig.adjrsquare;
        t = 1:0.1:12;
        hfit(cond) = plot(t,fitresult.a.*t.^fitresult.b + fitresult.c,'-k');
        hfit(cond).Color = h(cond).CData;%0.5*[1,1,1]
        hfit(cond).LineWidth = 2;

                    % Permute time and repeat the power law fit, to get a p-value for the
            % adjusted R^2 value
            rsquare_perms = nan(1,n_perms);
            adjrsquare_perms = nan(1,n_perms);
            for permI = 1:n_perms
                perm = randperm(12);
                y_perm = y(perm);
                [xData, yData] = prepareCurveData( x, y_perm);
                opts.Upper = [Inf 0 inf];
                [fitresult, gof] = fit( xData, yData, ft, opts ); 
                rsquare_perms(permI) = gof.rsquare;
                adjrsquare_perms(permI) = gof.adjrsquare;
            end
            prctile_num = sum(adjrsquare_perms<gof_orig.adjrsquare)/numel(adjrsquare_perms)*100;
            pval_str{cond}=['p=', num2str((100-prctile_num)/100)];
    end

end
ylims = ylim;
ylims(1) = floor(ylims(1)*10)/10;
ylims(2) = ceil(ylims(2)*10)/10;
ylim(ylims); % Round accuracy values to closest 10,20, etc.
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:12)
xlabel('session')
if plot_rt
    ylabel('RT [ms]')
    title('LT RT')
else
    ylabel('accuracy')
    title('Accuracy')
end

axis square
% legend([hfit],{['small R^2_a_d_j=' num2str(round(radj(1)*10)/10,2),' ',pval_str{1}],...
%     ['medium R^2_a_d_j='  num2str(round(radj(2)*10)/10,2),' ',pval_str{2}],...
%     ['large R^2_a_d_j=' num2str(round(radj(3)*10)/10,2),' ',pval_str{3}]})

set(gcf,'position',get(0,'screensize'))
if export_figures
    if plot_rt
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_LT_curves_RT [ms].png']);
    else
        fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_LT_curves_accuracy.png']);
    end
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end

    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end

%% Correlations between performance across tasks
rdm = mean(rdm_accuracy_per_week_and_coherence,3);
cd = mean(cd_accuracy_per_week_and_condition,3);
nc = nc_accuracy_per_week;
gng = mean(gng_accuracy_per_week_and_condition,3);
smb = mean(smb_reward_per_week_and_condition,3);
% Missing: LT, ITC

X = [rdm(:),cd(:),nc(:),gng(:),smb(:)];

correlations = corr(X,'rows','complete','type','Spearman');
for ii = 1:5
    correlations(ii,ii) = nan;
end
tasks = {'RDM','CD','NC','GNG','Bandit'}

figure('color','w')
imagesc(correlations)
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:5,'ytick',1:5,'xticklabel',tasks,'yticklabel',tasks)
title('Performance rank correlations')
subtitle('(calculated across all subjects)')
axis square
colorbar
colormap(coolwarm)
caxis([-0.4,0.4])

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_task_correlations',label,'.png']);
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end


%% Correlations between performance across tasks (calculated separately for each subject)
for sI = 1:90
    rdm = mean(rdm_accuracy_per_week_and_coherence(sI,:,:),3);
    cd = mean(cd_accuracy_per_week_and_condition(sI,:,:),3);
    nc = nc_accuracy_per_week(sI,:);
    gng = mean(gng_accuracy_per_week_and_condition(sI,:,:),3);
    smb = mean(smb_reward_per_week_and_condition(sI,:,:),3);
    X = [rdm(:),cd(:),nc(:),gng(:),smb(:)];
    correlations_per_subject(:,:,sI) = corr(X,'rows','complete','type','Spearman');
end

correlations_per_subject = fisherzinv(mean(fisherz(correlations_per_subject),3));

for ii = 1:5
    correlations_per_subject(ii,ii) = nan;
end
tasks = {'RDM','CD','NC','GNG','Bandit'}

figure('color','w')
nanimagesc(correlations_per_subject)
set(gca,'tickdir','out','fontsize',24,'box','off','xtick',1:5,'ytick',1:5,'xticklabel',tasks,'yticklabel',tasks)
title('Performance rank correlations')
subtitle('(average across individual correlations)')
axis square
colorbar
colormap(coolwarm)
caxis([-0.4,0.4])

set(gcf,'position',get(0,'screensize'))
if export_figures
    fig_file = fullfile(fig_dir,sub_dir,['hierarchical_completely_independent_task_correlations_within_subjects_',label,'.png']);
    if set_ylim
        fig_file = strrep(fig_file,'.png',['_ylim',num2str(ylims_all(1)),'_',num2str(ylims_all(2)),'.png']);
        ylim(ylims_all);
    end
    export_fig(gcf,fig_file,'-dpng','-r200');
    fig_file = strrep(fig_file,'.png','.pdf');
% fig_file = strrep(fig_file,'accuracy_rt_over_time','accuracy_rt_over_time/pdf');
exportgraphics(gcf,fig_file,'ContentType','vector')

end
