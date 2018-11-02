function get_cpm( im1_path, im2_path, flow_path, i_, j_, matchpath, model_path)

% setup 
dir_path = fileparts(mfilename('fullpath'));
piotr_path = fullfile(dir_path,'piotr_toolbox');
cpm_match_path = fullfile(dir_path,'build','cpm');

% addpath(genpath(SED_path));
addpath(genpath(piotr_path));

% input
% if (flow_path(1)~='/')
%     flow_path = fullfile(pwd,flow_path);
% end

% sparse correspondence estimation
match_path = fullfile(matchpath, sprintf('%s%s.txt', i_, j_));

cmd = sprintf('%s %s %s %s', cpm_match_path, im1_path, im2_path, match_path);
system(cmd);

% f = flow_path;

end % get_cpmflow

