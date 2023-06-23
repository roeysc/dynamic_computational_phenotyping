% % NOTE: THIS IS OLD CODE. I NOW USE:
% % /net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/shared_code/smb_precalculate_predictors/comppheno_get_smb_predictors.m




% addpath(genpath('/n/gershman_lab/users/Stan/1_code/roey_code/bandit'));
% clear
% 
% restart_kalman_with_correct_machine_assignment = false; % Restart the Kalman filter with the correct assignment of risky/safe variances (like in the original apar, where machines didn't change identity between blocks)
% restart_kalman_with_correct_machine_assignment_and_q = true; % This version also changes q (the observation variance), following https://github.com/sjgershm/exploration_uncertainty/blob/master/kalman_filter.m (Gershman, 2019)
% 
% data = compphenobandit_load_data();
% data_mat = data;
% data = data(:);
% idx_to_remove = [];
% for s = 1:length(data)
%     if data(s).N==0
%         idx_to_remove(end+1) = s;
%     end
% end
% data(idx_to_remove) = [];
% 
% %
% q0 = 36;           % observation variance
% count = 1;
% figure
% for q = [36]            % prior variance for initialization
%     % b = 100000; % Uncertainty bonus
%     % beta = 1; % Inverse temperature ?
%     % param = [q0,q,q,b,beta];
% %     param_risky0 = [q0,q,0]; % For blocks in which machine 0 is the risky machine
% %     param_risky1 = [q0,0,q]; % For blocks in which machine 1 is the risky machine
%     param_risky0 = [q,q,0.001]; % For blocks in which machine 0 is the risky machine THIS IS JUST TO SEE IF 1 GIVES BETTER RESULTS THAN 0 IN CASE 0 IS HARD FOR THE KALMAN FILTER
%     param_risky1 = [q0,0.001,q]; % For blocks in which machine 1 is the risky machine
%     %      param = [q0,0,q]; % Roey: Previously I had [q0,q,q] here but that's wrong, because the safe machine has 0 variance
%     param = [q0,q0/2,q0/2]; % Roey: This version implies that the subject aren't sure which of the two slot machines is the safe and which is the risky
%     
%     if restart_kalman_with_correct_machine_assignment
%         
%         for s = 1:length(data)
%             latents(s) = struct('m',[],'s',[]);
%             for b = 1:30
%                 idx = data(s).block==b;
%                 data_tmp = [];
%                 data_tmp.R = data(s).R(idx,:);
%                 data_tmp.block = data(s).block(idx);
%                 data_tmp.c = data(s).c(idx);
%                 data_tmp.r = data(s).r(idx);
%                 data_tmp.rt = data(s).rt(idx);
%                 data_tmp.trial = data(s).trial(idx);
%                 data_tmp.N = 10;
%                 data_tmp.C = 2;
%                 data_tmp.risky_machine = data(s).risky_machine(idx);
%                 if data_tmp.risky_machine(1)==0% &&& TODO: ROEY: THIS SHOULD BE 0, BUT I TRIED WITH 1 JUST TO SEE IF IT GIVES BETTER RESULTS
%                     latents_tmp = kalman_filter(param_risky0,data_tmp);
%                 else
%                     latents_tmp = kalman_filter(param_risky1,data_tmp);
%                 end
%                 latents(s).m = vertcat(latents(s).m,latents_tmp.m);
%                 latents(s).s = vertcat(latents(s).s,latents_tmp.s);
%             end
%         end
%         
%     elseif restart_kalman_with_correct_machine_assignment_and_q
%         
%             for s = 1:length(data)
%             latents(s) = struct('m',[],'s',[]);
%             for b = 1:30
%                 idx = data(s).block==b;
%                 data_tmp = [];
%                 data_tmp.R = data(s).R(idx,:);
%                 data_tmp.block = data(s).block(idx);
%                 data_tmp.c = data(s).c(idx);
%                 data_tmp.r = data(s).r(idx);
%                 data_tmp.rt = data(s).rt(idx);
%                 data_tmp.trial = data(s).trial(idx);
%                 data_tmp.N = 10;
%                 data_tmp.C = 2;
%                 data_tmp.risky_machine = data(s).risky_machine(idx);
%                 if data_tmp.risky_machine(1)==0% &&& TODO: ROEY: THIS SHOULD BE 0, BUT I TRIED WITH 1 JUST TO SEE IF IT GIVES BETTER RESULTS
%                     latents_tmp = kalman_filter2(param_risky0,data_tmp);
%                 else
%                     latents_tmp = kalman_filter2(param_risky1,data_tmp);
%                 end
%                 latents(s).m = vertcat(latents(s).m,latents_tmp.m);
%                 latents(s).s = vertcat(latents(s).s,latents_tmp.s);
%             end
%         end
%         
%     else % If we want to reset the initial values of the Kalman filter the same way regardless of the block type
%         for s = 1:length(data)
%             latents(s) = kalman_filter(param,data(s));
%         end
%     end
%     
%     clear b
%     X_all_subjects = [];
%     for s = 1:length(data) % loop over all sessions
%         V = latents(s).m(:,1) - latents(s).m(:,2); % V
%         RU = sqrt(latents(s).s(:,1)) - sqrt(latents(s).s(:,2)); % RU
%         TU = sqrt(latents(s).s(:,1) + latents(s).s(:,2)); % TU
%         C = double(data(s).c==1);
%         X = [V RU V./TU  ];
%         X_all_subjects = vertcat(X_all_subjects,X);
%         b(s,:) = glmfit(X,C,'binomial','link','probit','constant','off');
%     end
%     
%     %          subplot(1,6,count)
%     count = count+1;
%     %     idx = b(:,3)<100; % Remove one extreme values
%     idx = 1:size(b,1);
%     mu = mean(b(idx,:)); se = std(b(idx,:))./sqrt(size(b(idx,:),1));
%     barerrorbar(mu',se');
%     set(gca,'FontSize',25,'XTickLabel',{'V' 'RU' 'V/TU'});
%     title(['Prior var=' num2str(q)],'FontSize',25,'FontWeight','Bold');
%     ylabel('Coefficient','FontSize',25);
%     % ylim([-3 3])
%     axis square
% end
% 
% % Plot the input values for the probit regression
% binw = 1;
% figure
% histogram(X(:,1),'binw',binw)
% hold all
% histogram(X(:,2),'binw',binw)
% histogram(X(:,3),'binw',binw)
% legend({'V','RU','V/TU'})
% 
% 
% %% Plot the *weights* from the probit regression
% binw = 1;
% figure
% histogram(b(:,1),'binw',binw)
% hold all
% histogram(b(:,2),'binw',binw)
% histogram(b(:,3),'binw',binw)
% legend({'W_V','W_R_U','W_V_/_T_U'})
% 
% %% Compare with the values we get from Stan
% data_stan = readtable('/n/holyscratch01/gershman_lab/Lab/stan_output/smb/smb_mle_sessions_unbounded_3weights/combined_stan_outputs.csv');
% binw = 1;
% figure
% histogram(-data_stan.params_smb_1,'binw',binw)
% hold all
% histogram(-data_stan.params_smb_3,'binw',binw)
% histogram(-data_stan.params_smb_2,'binw',binw)
% legend({'W_V','W_R_U','W_V_/_T_U'})
% 
% 
% %%
% idx = 200;
% figure
% scatter(b(1:idx,1),-data_stan.params_smb_1(1:idx))
% identityLine(gca)
% 
% figure
% scatter(b(1:idx,2),-data_stan.params_smb_3(1:idx))
% identityLine(gca)
% 
% figure
% scatter(b(1:idx,3),-data_stan.params_smb_2(1:idx))
% identityLine(gca)
