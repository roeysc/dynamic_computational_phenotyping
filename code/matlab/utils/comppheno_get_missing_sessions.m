function missing_sessions = comppheno_get_missing_sessions(task)

%% Get missing sessions
% Get subjects in the order we always use them
subjects = comppheno_get_subjects();
data_dir = '/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/';

switch task
    case 'gng'
        fname = 'gng_data_for_stan_90s.csv';
    case 'rdm'
        fname = 'rdm_data_for_stan_90s.csv';
    case 'cd'
        fname = 'cd_data_for_stan_90s.csv';
    case 'nc'
        fname = 'nc_data_for_stan_90s.csv';
    case 'smb'
        fname = 'smb_data_for_stan_90s_updatedRegressorsNov22_rescaled_signVoverTU.csv';
    case 'lt'
        fname = 'lt_data_for_stan_90s.csv';
    case 'itc'
        fname = 'itc_data_for_stan_90s.csv';
end
        
% Load the tasl data
data = readtable(fullfile(data_dir,fname));

missing_sessions = zeros(length(subjects),12);
for sI = 1:length(subjects)
    sub = subjects{sI};
    for wI = 1:12
        idx = (strcmp(data.subjectId,sub) & data.weekId==wI);
        if all(idx==0)
            missing_sessions(sI,wI) = 1;
            continue
        end

    end
end
missing_sessions = logical(missing_sessions);