% tmp = readtable('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/4_hierarchical/gng_hierarchical_parallel/param_files/gng_parameters_powerlaw.csv')
% Run 
max_pc = 3;
out_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/exogenous_variables_test';
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

%% Exo PCA averaged across all 3 days of the week
session_to_hold_out = []; % Don't hold out any sessions
data_exo_pca = get_exo_inputs_after_PCA(session_to_hold_out,'');


csvwrite(fullfile(out_dir,'exo_PC1_valence.csv'),squeeze(data_exo_pca(1,:,:)))
csvwrite(fullfile(out_dir,'exo_PC2_arousal.csv'),squeeze(data_exo_pca(2,:,:)))
csvwrite(fullfile(out_dir,'exo_PC3_self_control.csv'),squeeze(data_exo_pca(3,:,:)))

%%
session_to_hold_out = []; % Don't hold out any sessions
data_exo_pca_day1 = get_exo_inputs_after_PCA(session_to_hold_out,'gng');
data_exo_pca_day1 = data_exo_pca_day1(1:max_pc,:,:);
data_exo_pca_day2 = get_exo_inputs_after_PCA(session_to_hold_out,'nc');
data_exo_pca_day2 = data_exo_pca_day2(1:max_pc,:,:);
data_exo_pca_day3 = get_exo_inputs_after_PCA(session_to_hold_out,'cd');
data_exo_pca_day3 = data_exo_pca_day3(1:max_pc,:,:);

%%
data_exo_pca = data_exo_pca_day1;
csvwrite(fullfile(out_dir,'exo_day1_gng_rdm_PC1_valence.csv'),squeeze(data_exo_pca(1,:,:)))
csvwrite(fullfile(out_dir,'exo_day1_gng_rdm_PC2_arousal.csv'),squeeze(data_exo_pca(2,:,:)))
csvwrite(fullfile(out_dir,'exo_day1_gng_rdm_PC3_self_control.csv'),squeeze(data_exo_pca(3,:,:)))

data_exo_pca = data_exo_pca_day2;
csvwrite(fullfile(out_dir,'exo_day2_lt_nc_smb_PC1_valence.csv'),squeeze(data_exo_pca(1,:,:)))
csvwrite(fullfile(out_dir,'exo_day2_lt_nc_smb_PC2_arousal.csv'),squeeze(data_exo_pca(2,:,:)))
csvwrite(fullfile(out_dir,'exo_day2_lt_nc_smb_PC3_self_control.csv'),squeeze(data_exo_pca(3,:,:)))

data_exo_pca = data_exo_pca_day3;
csvwrite(fullfile(out_dir,'exo_day3_cd_itc_PC1_valence.csv'),squeeze(data_exo_pca(1,:,:)))
csvwrite(fullfile(out_dir,'exo_day3_cd_itc_PC2_arousal.csv'),squeeze(data_exo_pca(2,:,:)))
csvwrite(fullfile(out_dir,'exo_day3_cd_itc_PC3_self_control.csv'),squeeze(data_exo_pca(3,:,:)))

%% Visualize individual subjects
pI = 1; 
tmp = squeeze(data_exo_pca_day1(pI,:,:));

figure('color','w')
num = 10; % Number of subjects to highlight
highlight = randperm(90);
highlight = highlight(1:num);
cmap = lines(num);
for sI = 1:90
    hold on
    if ismember(sI,highlight)
        continue
    else
        scatter(1:12,tmp(sI,:),'filled','MarkerFaceColor',[0,0,0],'MarkerFaceAlpha',0.1)
    end
end

count = 1;
for sI = 1:90
    hold on
    if ismember(sI,highlight)
        scatter(1:12,tmp(sI,:),'filled','MarkerFaceColor',cmap(count,:),'MarkerEdgeColor',[1,1,1])
        h = plot(1:12,tmp(sI,:));
h.Color = cmap(count,:);
        count = count+1;
    end
end

xlabel('session')
ylabel(['PC ' num2str(pI)])
set(gca,'tickdir','out')
