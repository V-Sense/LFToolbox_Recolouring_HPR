function spares_matches_centre(workingDir)

close all;

addpath(genpath('/utils/flow-code-matlab'));

MainDir = pwd;

MainImageSetDirPath = fullfile(MainDir,'/results', workingDir);

files = dir(MainImageSetDirPath);
dirFlags = [files.isdir];
images = files(~dirFlags);

FlowFileDirPath = fullfile(MainDir,'/flow', workingDir);

SparseMatchesDirPath = fullfile(MainDir, '/sparse_matches', workingDir);
mkdir(SparseMatchesDirPath);

tempPath = fullfile(FlowFileDirPath, '/forward_flow');
mkdir(tempPath);

% processing all non-black images - computing flow with centre image

for i = 1 : length(images)
    if( i == 1 || i == 2 || i == 3 || i == 13 || i == 14 || i == 15 || i == 16 || i == 30 ||... 			% black images
    	i == 196 || i == 210 || i == 211 || i == 212 || i == 213 || i == 223 || i == 224 || i == 225 ||... 	% black images
		i == 113 )																							% centre image
    	continue;
    end

    IMAGE1_PATH = fullfile(MainImageSetDirPath, '0707.png');
	IMAGE2_PATH = fullfile(MainImageSetDirPath, images(i).name);
    [~,IMAGE1_NAME,~] = fileparts(IMAGE1_PATH);
    [~,IMAGE2_NAME,~] = fileparts(IMAGE2_PATH);

	Forward_FLOW_SAVE_PATH = fullfile(FlowFileDirPath, 'forward_flow', [num2str(i-3),'.flo']);
	get_cpm(IMAGE1_PATH, IMAGE2_PATH, Forward_FLOW_SAVE_PATH, IMAGE1_NAME, IMAGE2_NAME, SparseMatchesDirPath);
end

end % function