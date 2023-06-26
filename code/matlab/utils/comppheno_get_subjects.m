function subjects  = comppheno_get_subjects(csv_file)
% Read the list of subjects from one of the behavioral data files and
% return it in the same order

comppheno_set_dirs % Load the comppheno_dir variable

if ~exist('csv_file','var')
    %csv_file = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/cd_data_for_stan_90s.csv'; % This file has the subjects in their correct order, and is therefor good to use.
    csv_file = fullfile(comppheno_data_dir,'cd_data_for_stan_90s.csv');
end
data = readtable(csv_file);
subjects = unique(data.subjectId,'stable');