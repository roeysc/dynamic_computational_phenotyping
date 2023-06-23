% comppheno_dir = 'SET/PATH/TO/THE/MAIN/MATLAB/DIRECTORYOF/THE/PROJECT'
comppheno_dir = '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/1_code/publish_code/code/matlab';
addpath(genpath(comppheno_dir));

comppheno_data_dir = fullfile(fileparts(fileparts(comppheno_dir)),'data');
comppheno_analysis_dir = fullfile(fileparts(fileparts(comppheno_dir)),'analysis');