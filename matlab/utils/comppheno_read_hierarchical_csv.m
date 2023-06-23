
function pheno_struct = comppheno_read_hierarchical_csv(csv_file,skip_priors)
% This read the csv file of different hierarchical models. The csv file is
% create by the python code OpenFile_roey.py

% csv_file = '/n/home10/rschurr/PycharmProjects/comp_phano_decision/gng_parameters_070822_gng88sub_rho20_hier.csv';
% task = 'gng';
if ~exist('skip_priors','var') || isempty(skip_priors)
    skip_priors = true;
end

% tmp = readtable(csv_file,'range','A46');
tmp = readtable(csv_file,'Delimiter',',');

parameters_list = {''}; % This will be filled automatically below

for rI = 2:size(tmp,1)
    name = strtrim(tmp{rI,1}); % strtrim makes sure to remove the new line that sometimes occures at the end of the file
    c = strsplit(name{1},'.');
%     if endsWith(c{1},'pr') %&& ~strcmp(c{1},'mu_pr')
%         continue
%     end
    if skip_priors && any(strfind(c{1},'_pr')) && ( ~any((strfind(c{1},'_pred')))) && ~any((strfind(c{1},'_prev'))) % This will ignore "prior" parameters (usually ending with _pr) but not predicted valued
        continue
    end
    if strcmp(c{1},'sigma') %&& ~strcmp(c{1},'mu_pr')
        continue
    end
    
    
% % %     if strcmp(c{1},'sigma') 
% % %         continue
% % %         eval([c{1} '(' c{2} ')=' num2str(tmp{rI,2}),';'])
% % %     end
% % %     if strcmp(c{1},'mu') || strcmp(c{1},'mu_pr') || strcmp(c{1},'mu_risk_lt') || strcmp(c{1},'mu_beta_lt')
% % %         continue
% % %         eval([c{1} '(' c{2} ')=' num2str(tmp{rI,2}),';'])
% % %     end
    
    learning_test = true;
    if learning_test
       if length(c)==1
            eval([c{1}, '=', num2str(tmp{rI,2}),';'])
       end
       if length(c)==2
            eval([c{1}, '(', c{2}, ')=', num2str(tmp{rI,2}),';'])
       end
       if length(c)==3
            eval([c{1} '(' c{2} ',' c{3} ')=' num2str(tmp{rI,2}),';'])
       end
       if length(c)==4
            eval([c{1} '(' c{2} ',' c{3} ,',', c{4}, ')=', num2str(tmp{rI,2}),';'])
       end
       parameters_list = unique([parameters_list(:);c{1}]);
        continue
    end
% % %     if ismember(task,{'cd'}) % These are tasks with 2 sub-indices (subejct x week x block)
% % %         eval([c{1} '(' c{2} ',' c{3} ',' c{4} ')=' num2str(tmp{rI,2}),';'])
% % %         parameters_list = unique([parameters_list(:);c{1}]);
% % %     elseif ismember(task,{'nc','itc','lt','gng','rdm'}) % These are tasks with 2 sub-indices (subejct x week)
% % % %     	if strfind(c{1},'delta_rdm_std') 
% % % %             continue % Just a silly fix for models where there was no delta_rdm_std parameter but I did include that as a title in the CSV file
% % % %         end
% % %         if strfind(c{1},'delta_rdm_noisy') 
% % %             continue % Just a silly fix for models where there was no delta_rdm_std parameter but I did include that as a title in the CSV file
% % %         end
% % %         if ~strfind(c{1},'delta_rdm_std') 
% % %             try
% % %                 eval([c{1} '(' c{2} ',' c{3} ')=' num2str(tmp{rI,2}),';'])
% % %             catch
% % %                 eval([c{1} '(' c{2} ',' c{3}(1:end-1) ')=' num2str(tmp{rI,2}),';']) % For some reason this CSV file has a weird symbol at the end of one cell, so this fixes it
% % %             end
% % %         else
% % %             eval([c{1} '(' c{2}  ')=' num2str(tmp{rI,2}),';'])
% % % 
% % %         end
% % %         parameters_list = unique([parameters_list(:);c{1}]);
% % %     elseif ismember(task,{'smb'}) % These are tasks with 2 sub-indices (subejct x week x paramNumber)
% % %         eval([c{1} '_' c{4} '(' c{2} ',' c{3} ')=' num2str(tmp{rI,2}),';'])
% % %         parameters_list = unique([parameters_list(:);{[c{1} '_' num2str(c{4})]}]);
% % %     end
end

parameters_list(cellfun(@isempty,parameters_list)) = []; % Remove the empty initial value

% Combine all phenotype matrices to a structure
for pI = 1:length(parameters_list)
    eval(['tmp = ', parameters_list{pI} ';'])
    pheno_struct.(parameters_list{pI}) = tmp;
end

