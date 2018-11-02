close all;
%% initialisation for colour transfer

addpath(genpath('data'));
addpath(genpath('results'));
addpath(genpath('sparse_matches'));

mg_initialiseMexFilesGT;
mg_initialiseMexFilesOMP;

%% initializing paths
%% PLACE HERE THE PATHS TO YOUR DATASETS

dataSet =   {
            % 'EPFL/Ankylosaurus_and_Diplodocus_1__Decoded';
            % 'EPFL/Bikes__Decoded';
            % 'EPFL/Color_Chart_1__Decoded';
            % 'EPFL/Danger_de_Mort__Decoded';
            % 'EPFL/Desktop__Decoded'
            % 'EPFL/Flowers__Decoded'
            % 'EPFL/Fountain_and_Vincent_2__Decoded'
            % 'EPFL/Friends_1__Decoded';
            % 'EPFL/ISO_Chart_12__Decoded';
            % 'EPFL/Magnets_1__Decoded';
            % 'EPFL/Stone_Pillars_Outside__Decoded';
            % 'EPFL/Vespa__Decoded';

            % 'INRIA/Bee_1__Decoded';
            % 'INRIA/Bee_2__Decoded';
            % 'INRIA/Bumblebee__Decoded';
            % 'INRIA/Cactus__Decoded';
            % 'INRIA/ChezEdgar__Decoded';
            % 'INRIA/Duck__Decoded';
            % 'INRIA/Field__Decoded';
            % 'INRIA/Fruits__Decoded';
            % 'INRIA/LensFlare__Decoded';
            % 'INRIA/Posts__Decoded';
            % 'INRIA/Rose__Decoded';
            % 'INRIA/SunnyRose__Decoded';
            % 'INRIA/Texture__Decoded';
            % 'INRIA/TinyMoon__Decoded';
            % 'INRIA/Translucent__Decoded';

            % 'Non_Lambertian/background__Decoded';
            % 'Non_Lambertian/bottle__Decoded';
            % 'Non_Lambertian/bracelet__Decoded';
            % 'Non_Lambertian/glasses1__Decoded';
            % 'Non_Lambertian/glasses2__Decoded';
            % 'Non_Lambertian/glasses3__Decoded';
            % 'Non_Lambertian/guiness__Decoded';
            % 'Non_Lambertian/Jam1__Decoded';
            % 'Non_Lambertian/Jam2__Decoded';
            % 'Non_Lambertian/Jam3__Decoded';
            % 'Non_Lambertian/teapot__Decoded';

            % 'V-SENSE/African_statue_2__Decoded';
            % 'V-SENSE/Cherry_tree_1__Decoded';
            % 'V-SENSE/Chicken_and_tea__Decoded';
            % 'V-SENSE/Oceane_outside_2__Decoded';
            % 'V-SENSE/Odette__Decoded';
            % 'V-SENSE/Raoul__Decoded';
            % 'V-SENSE/Rododendron_1__Decoded';
            % 'V-SENSE/Ukulele__Decoded';
            % 'V-SENSE/Wine_bottles_1__Decoded';
            
            };

for d = 1:length(dataSet)
    currentDir = dataSet{d, 1};
    matchDir = fullfile(pwd, '/sparse_matches', currentDir);
    matchFiles = dir(matchDir);
    dataDir = fullfile(pwd, '/data', currentDir);
    data = dir(dataDir);
    workingDir = fullfile(pwd, '/results', currentDir);
    mkdir(workingDir);
    mkdir(strcat(workingDir, '/param'));

    % performing hot pixel noise removal pre-recolouring
    hotPixelFix(30, 10, 7, dataDir, workingDir);

    % computing optical flow for the light fields (middle column + all rows, then between centre view and all views)
    sparse_matches_neighbours(currentDir); % modifications to be made in that file if specific processing necessary
    spares_matches_centre(currentDir); % same as above
    delete(fullfile(matchDir,'*.png')); % byproducts of the optical flow computation, no use for it

end

tStart = tic;
for d = 1:length(dataSet)
    currentDir = dataSet{d, 1};
    matchDir   = fullfile(pwd, '/sparse_matches', currentDir);
    matchFiles = dir(matchDir);
    dataDir    = fullfile(pwd, '/data', currentDir);
    data       = dir(dataDir);
    workingDir = fullfile(pwd, '/results', currentDir);
    files     = dir(workingDir); 
    dirFlags = [files.isdir];
    images = files(~dirFlags); % remove the folder names from the images list
    
    %% first section is for processing the middle column from centre to bottom

    for i = 1:length(images)
        if( ~isequal(images(i).name, '0707.png') && ~isequal(images(i).name, '0807.png') &&...
            ~isequal(images(i).name, '0907.png') && ~isequal(images(i).name, '1007.png') &&...
            ~isequal(images(i).name, '1107.png') && ~isequal(images(i).name, '1207.png') &&...
            ~isequal(images(i).name, '1307.png') )
            continue;  % ignoring files '.', '..', and black images (always in the same positions thank God)
        end

        % parsing images (palette + target)
        
        cpalette_file       = fullfile(workingDir, '0707.png');
        npalette_file       = fullfile(workingDir, images(i).name);
        target_file         = fullfile(workingDir, images(i+15).name);
        param_file_new      = fullfile(workingDir, 'param', strcat(images(i+15).name(1:end-4), '.mat'));
        param_file_pal      = fullfile(workingDir, 'param', strcat(images(i).name(1:end-4), '.mat'));
        [~,cpalette_name,~] = fileparts(cpalette_file);
        [~,npalette_name,~] = fileparts(npalette_file);
        [~,target_name,~]   = fileparts(target_file);
        cpalette            = imread(cpalette_file);
        npalette            = imread(npalette_file);
        target              = imread(target_file);
        
        fprintf('Current folder : %s.\n', currentDir);
        fprintf('Processing palette : %s.png and target : %s.png with centre image.\n', npalette_name, target_name);

        % using the matches to find the corresponding colours (in CIELab space) of the pixels in the LF images
        cmatch_file = sprintf('%s%s.txt', cpalette_name, target_name);
        cmatch_path = fullfile(matchDir, cmatch_file);
        cmatch      = dlmread(cmatch_path);

        nmatch_file = sprintf('%s%s.txt', npalette_name, target_name);
        nmatch_path = fullfile(matchDir, nmatch_file);
        nmatch      = dlmread(nmatch_path);
        
        cmatch_size = size(cmatch, 1);
        nmatch_size = size(nmatch, 1);
        
        colours_palette(1:cmatch_size+nmatch_size, 3) = double(0); % matched colours in image 1 (used by CT)
        colours_target(1:cmatch_size+nmatch_size, 3)  = double(0); % matched colours in image 2

        for k = 1:cmatch_size
            colours_palette(k,:) = cpalette(cmatch(k,2) + 1, cmatch(k,1) + 1, :);
            colours_target(k,:)  = target(cmatch(k,4) + 1, cmatch(k,3) + 1, :);
        end

        for k = cmatch_size+1:cmatch_size+nmatch_size
            colours_palette(k,:) = npalette(nmatch(k-cmatch_size,2) + 1, nmatch(k-cmatch_size,1) + 1, :);
            colours_target(k,:)  = target(nmatch(k-cmatch_size,4) + 1, nmatch(k-cmatch_size,3) + 1, :);
        end

        % using these colours, applying colour transfer between palette image and target image
        %load previous parameters
        if(exist(param_file_pal))
            load(param_file_pal);
        else 
            param = [];
        end
        
        tic
        [param] = ctfunction_corr(cpalette, target, 8, 2, colours_palette, colours_target, target_file, param, 'CIELab');
        toc
        save(param_file_new, 'param');
    end

    %% second section is for processing the middle column from centre to top

    for i = length(images):-1:1
        if( ~isequal(images(i).name, '0707.png') && ~isequal(images(i).name, '0607.png') &&...
            ~isequal(images(i).name, '0507.png') && ~isequal(images(i).name, '0407.png') &&...
            ~isequal(images(i).name, '0307.png') && ~isequal(images(i).name, '0207.png') &&...
            ~isequal(images(i).name, '0107.png') )
            continue;
        end
    
        cpalette_file       = fullfile(workingDir, '0707.png');
        npalette_file       = fullfile(workingDir, images(i).name);
        target_file         = fullfile(workingDir, images(i-15).name);
        param_file_new      = fullfile(workingDir, 'param', strcat(images(i-15).name(1:end-4), '.mat'));
        param_file_pal      = fullfile(workingDir, 'param', strcat(images(i).name(1:end-4), '.mat'));
        [~,cpalette_name,~] = fileparts(cpalette_file);
        [~,npalette_name,~] = fileparts(npalette_file);
        [~,target_name,~]   = fileparts(target_file);
        cpalette            = imread(cpalette_file);
        npalette            = imread(npalette_file);
        target              = imread(target_file);

        fprintf('Current folder : %s.\n', currentDir);
        fprintf('Processing palette : %s.png and target : %s.png with centre image.\n', npalette_name, target_name);

        cmatch_file = sprintf('%s%s.txt', cpalette_name, target_name);
        cmatch_path = fullfile(matchDir, cmatch_file);
        cmatch      = dlmread(cmatch_path);

        nmatch_file = sprintf('%s%s.txt', target_name, npalette_name);
        nmatch_path = fullfile(matchDir, nmatch_file);
        nmatch      = dlmread(nmatch_path);
        
        cmatch_size = size(cmatch, 1);
        nmatch_size = size(nmatch, 1);

        colours_palette(1:cmatch_size+nmatch_size, 3) = double(0);
        colours_target(1:cmatch_size+nmatch_size, 3) = double(0);

        for k = 1:cmatch_size
            colours_palette(k,:) = cpalette(cmatch(k,2) + 1, cmatch(k,1) + 1, :);
            colours_target(k,:)  = target(cmatch(k,4) + 1, cmatch(k,3) + 1, :);
        end

        for k = cmatch_size+1:cmatch_size+nmatch_size
            colours_palette(k,:) = npalette(nmatch(k-cmatch_size,4) + 1, nmatch(k-cmatch_size,3) + 1, :);
            colours_target(k,:)  = target(nmatch(k-cmatch_size,2) + 1, nmatch(k-cmatch_size,1) + 1, :);
        end

        if(exist(param_file_pal))
            load(param_file_pal);
        else 
            param = [];
        end
        
        tic
        [param] = ctfunction_corr(cpalette, target, 8, 2, colours_palette, colours_target, target_file, param, 'CIELab');
        toc
        save(param_file_new, 'param');
    end

    %% third section for processing a first time the dark images in the corners (not the black ones) on the right

    for i = 1:length(images)
        if( ~isequal(images(i).name, '0010.png') && ~isequal(images(i).name, '0112.png') &&...
            ~isequal(images(i).name, '0213.png') && ~isequal(images(i).name, '0313.png') &&...    
            ~isequal(images(i).name, '1113.png') && ~isequal(images(i).name, '1213.png') &&...    
            ~isequal(images(i).name, '1312.png') && ~isequal(images(i).name, '1410.png') )
            continue;
        end
            
        palette_file        = fullfile(workingDir, images(i).name);
        target_file         = fullfile(workingDir, images(i+1).name);
        [~,palette_name,~]  = fileparts(palette_file);
        [~,target_name,~]   = fileparts(target_file);
        palette             = imread(palette_file);
        target              = imread(target_file);
        
        fprintf('Current folder : %s.\n', currentDir);
        fprintf('Processing palette : %s.png and target : %s.png.\n', palette_name, target_name);

        match_file = sprintf('%s%s.txt', palette_name, target_name);
        match_path = fullfile(matchDir, match_file);
        match      = dlmread(match_path);
        
        match_size = size(match, 1);

        colours_palette(1:match_size, 3) = double(0);
        colours_target(1:match_size, 3)  = double(0);

        for k = 1:match_size
            colours_palette(k,:) = palette(match(k,2) + 1, match(k,1) + 1, :);
            colours_target(k,:)  = target(match(k,4) + 1, match(k,3) + 1, :);
        end

        tic
        ctfunction_corr(palette, target, 8, 2, colours_palette, colours_target, target_file,[], 'RGB');
        toc
    end

    %% fourth section for processing a first time the dark images in the corners (not the black ones) on the left

    for i = length(images):-1:1
        if( ~isequal(images(i).name, '0004.png') && ~isequal(images(i).name, '0102.png') &&...
            ~isequal(images(i).name, '0201.png') && ~isequal(images(i).name, '0301.png') &&...    
            ~isequal(images(i).name, '1101.png') && ~isequal(images(i).name, '1201.png') &&...    
            ~isequal(images(i).name, '1302.png') && ~isequal(images(i).name, '1404.png') )
            continue;
        end             
        
        palette_file       = fullfile(workingDir, images(i).name);
        target_file        = fullfile(workingDir, images(i-1).name);
        [~,palette_name,~] = fileparts(palette_file);
        [~,target_name,~]  = fileparts(target_file);
        palette            = imread(palette_file);
        target             = imread(target_file);
        
        fprintf('Current folder : %s.\n', currentDir);
        fprintf('Processing palette : %s.png and target : %s.png.\n', palette_name, target_name);

        match_file = sprintf('%s%s.txt', target_name, palette_name);
        match_path = fullfile(matchDir, match_file);
        match      = dlmread(match_path);
        
        match_size = size(match, 1);

        colours_palette(1:match_size, 3) = double(0);
        colours_target(1:match_size, 3)  = double(0);

        for k = 1:match_size
            colours_palette(k,:) = palette(match(k,2) + 1, match(k,1) + 1, :);
            colours_target(k,:)  = target(match(k,4) + 1, match(k,3) + 1, :);
        end

        tic
        ctfunction_corr(palette, target, 8, 2, colours_palette, colours_target, target_file, [], 'RGB');
        toc
    end

    %% fifth section to process rows from centre to right

    for i = 1:length(images)
        if( ~isequal(images(i).name, '0007.png') && ~isequal(images(i).name, '0008.png') &&...
            ~isequal(images(i).name, '0009.png') && ~isequal(images(i).name, '0010.png') &&...    
            ~isequal(images(i).name, '0107.png') && ~isequal(images(i).name, '0108.png') &&...
            ~isequal(images(i).name, '0109.png') && ~isequal(images(i).name, '0110.png') &&...
            ~isequal(images(i).name, '0111.png') && ~isequal(images(i).name, '0112.png') &&...
            ~isequal(images(i).name, '0207.png') && ~isequal(images(i).name, '0208.png') &&...
            ~isequal(images(i).name, '0209.png') && ~isequal(images(i).name, '0210.png') &&...
            ~isequal(images(i).name, '0211.png') && ~isequal(images(i).name, '0212.png') &&...
            ~isequal(images(i).name, '0213.png') && ~isequal(images(i).name, '0307.png') &&...
            ~isequal(images(i).name, '0308.png') && ~isequal(images(i).name, '0309.png') &&...
            ~isequal(images(i).name, '0310.png') && ~isequal(images(i).name, '0311.png') &&...
            ~isequal(images(i).name, '0312.png') && ~isequal(images(i).name, '0313.png') &&...
            ~isequal(images(i).name, '0407.png') && ~isequal(images(i).name, '0408.png') &&...
            ~isequal(images(i).name, '0409.png') && ~isequal(images(i).name, '0410.png') &&...
            ~isequal(images(i).name, '0411.png') && ~isequal(images(i).name, '0412.png') &&...
            ~isequal(images(i).name, '0413.png') && ~isequal(images(i).name, '0507.png') &&...
            ~isequal(images(i).name, '0508.png') && ~isequal(images(i).name, '0509.png') &&...
            ~isequal(images(i).name, '0510.png') && ~isequal(images(i).name, '0511.png') &&...
            ~isequal(images(i).name, '0512.png') && ~isequal(images(i).name, '0513.png') &&...
            ~isequal(images(i).name, '0607.png') && ~isequal(images(i).name, '0608.png') &&...
            ~isequal(images(i).name, '0609.png') && ~isequal(images(i).name, '0610.png') &&...
            ~isequal(images(i).name, '0611.png') && ~isequal(images(i).name, '0612.png') &&...
            ~isequal(images(i).name, '0613.png') && ~isequal(images(i).name, '0707.png') &&...
            ~isequal(images(i).name, '0708.png') && ~isequal(images(i).name, '0709.png') &&...
            ~isequal(images(i).name, '0710.png') && ~isequal(images(i).name, '0711.png') &&...
            ~isequal(images(i).name, '0712.png') && ~isequal(images(i).name, '0713.png') &&...
            ~isequal(images(i).name, '0807.png') && ~isequal(images(i).name, '0808.png') &&...
            ~isequal(images(i).name, '0809.png') && ~isequal(images(i).name, '0810.png') &&...
            ~isequal(images(i).name, '0811.png') && ~isequal(images(i).name, '0812.png') &&...
            ~isequal(images(i).name, '0813.png') && ~isequal(images(i).name, '0907.png') &&...
            ~isequal(images(i).name, '0908.png') && ~isequal(images(i).name, '0909.png') &&...
            ~isequal(images(i).name, '0910.png') && ~isequal(images(i).name, '0911.png') &&...
            ~isequal(images(i).name, '0912.png') && ~isequal(images(i).name, '0913.png') &&...
            ~isequal(images(i).name, '1007.png') && ~isequal(images(i).name, '1008.png') &&...
            ~isequal(images(i).name, '1009.png') && ~isequal(images(i).name, '1010.png') &&...
            ~isequal(images(i).name, '1011.png') && ~isequal(images(i).name, '1012.png') &&...
            ~isequal(images(i).name, '1013.png') && ~isequal(images(i).name, '1107.png') &&...
            ~isequal(images(i).name, '1108.png') && ~isequal(images(i).name, '1109.png') &&...
            ~isequal(images(i).name, '1110.png') && ~isequal(images(i).name, '1111.png') &&...
            ~isequal(images(i).name, '1112.png') && ~isequal(images(i).name, '1113.png') &&...
            ~isequal(images(i).name, '1207.png') && ~isequal(images(i).name, '1208.png') &&...
            ~isequal(images(i).name, '1209.png') && ~isequal(images(i).name, '1210.png') &&...
            ~isequal(images(i).name, '1211.png') && ~isequal(images(i).name, '1212.png') &&...
            ~isequal(images(i).name, '1213.png') && ~isequal(images(i).name, '1307.png') &&...
            ~isequal(images(i).name, '1308.png') && ~isequal(images(i).name, '1309.png') &&...
            ~isequal(images(i).name, '1310.png') && ~isequal(images(i).name, '1311.png') &&...
            ~isequal(images(i).name, '1312.png') && ~isequal(images(i).name, '1407.png') &&...
            ~isequal(images(i).name, '1408.png') && ~isequal(images(i).name, '1409.png') &&...
            ~isequal(images(i).name, '1410.png') )
            continue;
        end

        cpalette_file       = fullfile(workingDir, '0707.png');
        npalette_file       = fullfile(workingDir, images(i).name);
        target_file         = fullfile(workingDir, images(i+1).name);
        param_file_new      = fullfile(workingDir, 'param', strcat(images(i+1).name(1:end-4), '.mat'));
        param_file_pal      = fullfile(workingDir, 'param', strcat(images(i).name(1:end-4), '.mat'));
        [~,cpalette_name,~] = fileparts(cpalette_file);
        [~,npalette_name,~] = fileparts(npalette_file);
        [~,target_name,~]   = fileparts(target_file);
        cpalette            = imread(cpalette_file);
        npalette            = imread(npalette_file);
        target              = imread(target_file);
        
        fprintf('Current folder : %s.\n', currentDir);
        fprintf('Processing palette : %s.png and target : %s.png with centre image.\n', npalette_name, target_name);

        cmatch_file = sprintf('%s%s.txt', cpalette_name, target_name);
        cmatch_path = fullfile(matchDir, cmatch_file);
        cmatch      = dlmread(cmatch_path);

        nmatch_file = sprintf('%s%s.txt', npalette_name, target_name);
        nmatch_path = fullfile(matchDir, nmatch_file);
        nmatch      = dlmread(nmatch_path);
        
        cmatch_size = size(cmatch, 1);
        nmatch_size = size(nmatch, 1);

        colours_palette(1:cmatch_size+nmatch_size, 3) = double(0);
        colours_target(1:cmatch_size+nmatch_size, 3)  = double(0);

        for k = 1:cmatch_size
            colours_palette(k,:) = cpalette(cmatch(k,2) + 1, cmatch(k,1) + 1, :);
            colours_target(k,:)  = target(cmatch(k,4) + 1, cmatch(k,3) + 1, :);
        end

        for k = cmatch_size+1:cmatch_size+nmatch_size
            colours_palette(k,:) = npalette(nmatch(k-cmatch_size,2) + 1, nmatch(k-cmatch_size,1) + 1, :);
            colours_target(k,:)  = target(nmatch(k-cmatch_size,4) + 1, nmatch(k-cmatch_size,3) + 1, :);
        end
        
        if(exist(param_file_pal))
            load(param_file_pal);
        else 
            param = [];
        end
        
        tic
        [param] = ctfunction_corr(cpalette, target, 8, 1, colours_palette, colours_target, target_file, param, 'CIELab');
        toc
        save(param_file_new, 'param');
    end

    %% sixth section to process rows from centre to left

    for i = length(images):-1:1
        if( ~isequal(images(i).name, '0007.png') && ~isequal(images(i).name, '0006.png') &&...
            ~isequal(images(i).name, '0005.png') && ~isequal(images(i).name, '0004.png') &&...
            ~isequal(images(i).name, '0107.png') && ~isequal(images(i).name, '0106.png') &&...
            ~isequal(images(i).name, '0105.png') && ~isequal(images(i).name, '0104.png') &&...
            ~isequal(images(i).name, '0103.png') && ~isequal(images(i).name, '0102.png') &&...
            ~isequal(images(i).name, '0207.png') && ~isequal(images(i).name, '0206.png') &&...
            ~isequal(images(i).name, '0205.png') && ~isequal(images(i).name, '0204.png') &&...
            ~isequal(images(i).name, '0203.png') && ~isequal(images(i).name, '0202.png') &&...
            ~isequal(images(i).name, '0201.png') && ~isequal(images(i).name, '0307.png') &&...
            ~isequal(images(i).name, '0306.png') && ~isequal(images(i).name, '0305.png') &&...
            ~isequal(images(i).name, '0304.png') && ~isequal(images(i).name, '0303.png') &&...
            ~isequal(images(i).name, '0302.png') && ~isequal(images(i).name, '0301.png') &&...
            ~isequal(images(i).name, '0407.png') && ~isequal(images(i).name, '0406.png') &&...
            ~isequal(images(i).name, '0405.png') && ~isequal(images(i).name, '0404.png') &&...
            ~isequal(images(i).name, '0403.png') && ~isequal(images(i).name, '0402.png') &&...
            ~isequal(images(i).name, '0401.png') && ~isequal(images(i).name, '0507.png') &&...
            ~isequal(images(i).name, '0506.png') && ~isequal(images(i).name, '0505.png') &&...
            ~isequal(images(i).name, '0504.png') && ~isequal(images(i).name, '0503.png') &&...
            ~isequal(images(i).name, '0502.png') && ~isequal(images(i).name, '0501.png') &&...
            ~isequal(images(i).name, '0607.png') && ~isequal(images(i).name, '0606.png') &&...
            ~isequal(images(i).name, '0605.png') && ~isequal(images(i).name, '0604.png') &&...
            ~isequal(images(i).name, '0603.png') && ~isequal(images(i).name, '0602.png') &&...
            ~isequal(images(i).name, '0601.png') && ~isequal(images(i).name, '0707.png') &&...
            ~isequal(images(i).name, '0706.png') && ~isequal(images(i).name, '0705.png') &&...
            ~isequal(images(i).name, '0704.png') && ~isequal(images(i).name, '0703.png') &&...
            ~isequal(images(i).name, '0702.png') && ~isequal(images(i).name, '0701.png') &&...
            ~isequal(images(i).name, '0807.png') && ~isequal(images(i).name, '0806.png') &&...
            ~isequal(images(i).name, '0805.png') && ~isequal(images(i).name, '0804.png') &&...
            ~isequal(images(i).name, '0803.png') && ~isequal(images(i).name, '0802.png') &&...
            ~isequal(images(i).name, '0801.png') && ~isequal(images(i).name, '0907.png') &&...
            ~isequal(images(i).name, '0906.png') && ~isequal(images(i).name, '0905.png') &&...
            ~isequal(images(i).name, '0904.png') && ~isequal(images(i).name, '0903.png') &&...
            ~isequal(images(i).name, '0902.png') && ~isequal(images(i).name, '0901.png') &&...
            ~isequal(images(i).name, '1007.png') && ~isequal(images(i).name, '1006.png') &&...
            ~isequal(images(i).name, '1005.png') && ~isequal(images(i).name, '1004.png') &&...
            ~isequal(images(i).name, '1003.png') && ~isequal(images(i).name, '1002.png') &&...
            ~isequal(images(i).name, '1001.png') && ~isequal(images(i).name, '1107.png') &&...
            ~isequal(images(i).name, '1106.png') && ~isequal(images(i).name, '1105.png') &&...
            ~isequal(images(i).name, '1104.png') && ~isequal(images(i).name, '1103.png') &&...
            ~isequal(images(i).name, '1102.png') && ~isequal(images(i).name, '1101.png') &&...
            ~isequal(images(i).name, '1207.png') && ~isequal(images(i).name, '1206.png') &&...
            ~isequal(images(i).name, '1205.png') && ~isequal(images(i).name, '1204.png') &&...
            ~isequal(images(i).name, '1203.png') && ~isequal(images(i).name, '1202.png') &&...
            ~isequal(images(i).name, '1201.png') && ~isequal(images(i).name, '1307.png') &&...
            ~isequal(images(i).name, '1306.png') && ~isequal(images(i).name, '1305.png') &&...
            ~isequal(images(i).name, '1304.png') && ~isequal(images(i).name, '1303.png') &&...
            ~isequal(images(i).name, '1302.png') && ~isequal(images(i).name, '1407.png') &&...
            ~isequal(images(i).name, '1406.png') && ~isequal(images(i).name, '1405.png') &&...
            ~isequal(images(i).name, '1404.png') )
            continue;
        end

        cpalette_file       = fullfile(workingDir, '0707.png');
        npalette_file       = fullfile(workingDir, images(i).name);
        target_file         = fullfile(workingDir, images(i-1).name);
        param_file_new      = fullfile(workingDir, 'param', strcat(images(i-1).name(1:end-4), '.mat'));
        param_file_pal      = fullfile(workingDir, 'param', strcat(images(i).name(1:end-4), '.mat'));
        [~,cpalette_name,~] = fileparts(cpalette_file);
        [~,npalette_name,~] = fileparts(npalette_file);
        [~,target_name,~]   = fileparts(target_file);
        cpalette            = imread(cpalette_file);
        npalette            = imread(npalette_file);
        target              = imread(target_file);

        fprintf('Current folder : %s.\n', currentDir);
        fprintf('Processing palette : %s.png and target : %s.png with centre image.\n', npalette_name, target_name);

        cmatch_file = sprintf('%s%s.txt', cpalette_name, target_name);
        cmatch_path = fullfile(matchDir, cmatch_file);
        cmatch      = dlmread(cmatch_path);

        nmatch_file = sprintf('%s%s.txt', target_name, npalette_name);
        nmatch_path = fullfile(matchDir, nmatch_file);
        nmatch      = dlmread(nmatch_path);
        
        cmatch_size = size(cmatch, 1);
        nmatch_size = size(nmatch, 1);

        colours_palette(1:cmatch_size+nmatch_size, 3) = double(0);
        colours_target(1:cmatch_size+nmatch_size, 3)  = double(0);

        for k = 1:cmatch_size
            colours_palette(k,:) = cpalette(cmatch(k,2) + 1, cmatch(k,1) + 1, :);
            colours_target(k,:)  = target(cmatch(k,4) + 1, cmatch(k,3) + 1, :);
        end

        for k = cmatch_size+1:cmatch_size+nmatch_size
            colours_palette(k,:) = npalette(nmatch(k-cmatch_size,4) + 1, nmatch(k-cmatch_size,3) + 1, :);
            colours_target(k,:)  = target(nmatch(k-cmatch_size,2) + 1, nmatch(k-cmatch_size,1) + 1, :);
        end
        
        if(exist(param_file_pal))
            load(param_file_pal);
        else 
            param = [];
        end
        
        tic
        [param] = ctfunction_corr(cpalette, target, 8, 2, colours_palette, colours_target, target_file, param, 'CIELab');
        toc
        save(param_file_new, 'param');
    end

    %% some minor cleanup

    delete(fullfile(pwd, '/tmp/', '*'));

    % performing hot pixel noise removal post-recolouring for dark corner images
    hotPixelFix_corners(30, 10, 7, workingDir);
    
end % for loop (on folders)

tElapsed = toc(tStart)
fprintf('Total processing time = %f.\n', tElapsed);