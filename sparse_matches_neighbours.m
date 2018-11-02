function sparse_matches_neighbours(workingDir)

close all;

addpath(genpath('/utils/flow-code-matlab'));

MainDir = pwd;

MainImageSetDirPath = fullfile(MainDir,'/results', workingDir); % same as working dir in main code

files = dir(MainImageSetDirPath);
dirFlags = [files.isdir];
images = files(~dirFlags);

FlowFileDirPath = fullfile(MainDir,'/flow', workingDir);

SparseMatchesDirPath = fullfile(MainDir, '/sparse_matches', workingDir);
mkdir(SparseMatchesDirPath);

tempPath = fullfile(FlowFileDirPath, '/forward_flow');
mkdir(tempPath);

% processing all rows

for i = 1 : length(images)
    if( i == 1 || i == 2 || i == 3 || i == 13 || i == 14 || i == 15 || i == 16 || i == 30 ||... 			% black images
    	i == 196 || i == 210 || i == 211 || i == 212 || i == 213 || i == 223 || i == 224 || i == 225 ||... 	% black images
    	i+1 == 13 || i+1 == 30 || i+1 == 210 || i+1 == 223 ||... 											% following image is black
    	i == 45 || i == 60 || i == 75 || i == 90 || i == 105 || i == 120 ||... 								% end of row
    	i == 135 || i == 150 || i == 165 || i == 180 || i == 195 )											% end of row
    	continue;
    end

	IMAGE1_PATH = fullfile(MainImageSetDirPath, images(i).name);
    IMAGE2_PATH = fullfile(MainImageSetDirPath, images(i+1).name);
    [~,IMAGE1_NAME,~] = fileparts(IMAGE1_PATH);
    [~,IMAGE2_NAME,~] = fileparts(IMAGE2_PATH);

	Forward_FLOW_SAVE_PATH = fullfile(FlowFileDirPath, 'forward_flow', [num2str(i-3),'.flo']);
	get_cpm(IMAGE1_PATH, IMAGE2_PATH, Forward_FLOW_SAVE_PATH, IMAGE1_NAME, IMAGE2_NAME, SparseMatchesDirPath);
end

% processing centre column

for i = 1 : length(images)
    if( i ~= 8 && i ~= 23 && i ~= 38 && i ~= 53 && i ~= 68 && i ~= 83 && i ~= 98 &&...
    	i ~= 113 && i ~= 128 && i ~= 143 && i ~= 158 && i ~= 173 && i ~= 188 && i ~= 203 )
    	continue;
    end
    
    IMAGE1_PATH = fullfile(MainImageSetDirPath, images(i).name);
    IMAGE2_PATH = fullfile(MainImageSetDirPath, images(i+15).name);
    [~,IMAGE1_NAME,~] = fileparts(IMAGE1_PATH);
    [~,IMAGE2_NAME,~] = fileparts(IMAGE2_PATH);

    Forward_FLOW_SAVE_PATH = fullfile(FlowFileDirPath, 'forward_flow', [num2str(i-3),'col.flo']);
    get_cpm(IMAGE1_PATH, IMAGE2_PATH, Forward_FLOW_SAVE_PATH, IMAGE1_NAME, IMAGE2_NAME, SparseMatchesDirPath);
end

end % function