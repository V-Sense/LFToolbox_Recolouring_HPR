function [] = ctfunction_parallel(paletteFile, targetFolder, newPalette, bitCoding, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ctfunction(paletteFile, targetFile, clusterFun, nColors,colourSpace )
% Computes the colour transfer result between a target image and palette
% image that are aligned, eg. pixels at the same location in both images
% are correspondences.
%
% paletteFile:  The name of the palette file to be processed. The colour
%               distribution of the palette file will be mapped to the target image. 
% targetFile:   The name of the target file to be processed.The colour
%               distribution of the palette file will be mapped to the target image.
% clusterFun:   'KMeans' or 'MVQ'(default). The clustering function used.
% nColors:      The number of clusters computed by the clustering
%               function. The default is 50.
% colourSpace:  The colour space in which the registration is performed. The default is 'RGB'. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if(nargin < 4)
    %Four input arguements must be given
    error('Four input arguements are required: palette image, target folder, first target image, encoding of target images');
elseif(nargin == 4)
    %If a clustering function is not specified, use Matlab's MVQ fucntion.
    %If the number of clusters to be used is not specified, use 50. 
    clusterFun = 'MVQ';
    nColors = 50;
    colourSpace = 'RGB';
elseif(nargin == 5)
    clusterFun = varargin{1};
    nColors = 50;
    colourSpace = 'RGB';
elseif(nargin == 6)
    clusterFun = varargin{1};
    nColors = varargin{2};
    colourSpace = 'RGB';
elseif(nargin == 7)
    clusterFun = varargin{1};
    nColors = varargin{2};
    colourSpace = varargin{3};
end
addpath(genpath('../L2RegistrationForCT'));
    
close all;

%read target image and reshape
LFfolder = targetFolder;

if ~isdir(LFfolder)
    errorMessage = sprintf('Error : following folder does not exist in this realm : \n%s', LFfolder);
    uiwait(warndlg(errorMessage));
    return;
end

filePattern = fullfile(LFfolder, '*.png');
lfield = dir(filePattern);

he1 = imread(newPalette);

%imshow(he1);

% if(bitCoding == 8)
%     he1 = im2uint8(he1);
% end

if(strcmp(colourSpace, 'CIELab') && strcmp(clusterFun, 'KMeans'))
	he1	= rgb2lab(he1);
elseif(strcmp(colourSpace, 'CIELab') && strcmp(clusterFun, 'MVQ'))
    disp('MVQ clustering cannot be applied in CIELAB space in this implementation. Using RGB space instead.');
    colourSpace = 'RGB';
end

dfull1 = double(he1);
fnrows = size(dfull1,1);
fncols = size(dfull1,2);
full_transform = reshape(dfull1,fnrows*fncols,3);

%read palette image
he2 = imread(paletteFile);

%imshow(he2);

if(bitCoding == 16)
    he2 = im2uint16(he2);
end
if(strcmp(colourSpace, 'CIELab'))
	he2	= rgb2lab(he2);
end

%cluster target and palette image to get most dominant colours
disp('Clustering of palette and target image started...');
switch(clusterFun)
    case 'KMeans'
        X = mg_applyKMeans(he1,nColors);
        Y = mg_applyKMeans(he2,nColors);
    case 'MVQ'
        if(bitCoding == 8)
            X = mg_quantImage8(he1, nColors);
            Y = mg_quantImage8(he2, nColors);
        elseif(bitCoding == 16)
            X = mg_quantImage16(he1, nColors);
            Y = mg_quantImage16(he2, nColors);
        else
            disp('Image is either 8 or 16-bit, you fool!');
        end
end
disp('Clustering of palette and target image finished.');

%initialise some parameters used to control the registration 
[config] = mg_initialize_config(X,Y, colourSpace, bitCoding);
%save 'configbef.mat' 'config';
disp('Registration of colours started...');
%register the colours using an annealling scheme, estimate a TPS transformation 
for i = 1:config.AnnSteps
    config.iter = (config.AnnSteps-i+1);
    [param, ~, ~, config] = gmmreg_rbf_L2(config);
     if(i ~= config.AnnSteps)
        config.scale = .5*config.scale;
    end
end
disp('Registration of colours finished.');
%save 'configaf.mat' 'config';

%saving parameters on file for subsequent use in the loop (hopefully prevents errors?)
%save('config.mat', 'config');
%save('param.mat', 'param');

%apply the colour transformation to the target image to get the result
fprintf(1, 'Applying colour transfer to image : %s\n', newPalette);
fullTransform = mg_transform_tps_parallel(param, full_transform, config.ctrl_pts);
disp('Finished.');
if(strcmp(colourSpace, 'CIELab'))
       fullTransform = lab2rgb(fullTransform);
       ind1 = find(fullTransform < 0);
       fullTransform(ind1) = 0;
       ind2 = find(fullTransform > 1);
       fullTransform(ind2) = 1;
       fullTransform = 255*fullTransform;
end

finalResult = reshape(fullTransform, [fnrows fncols 3]);
if(bitCoding == 8)
    finalResult = uint8(finalResult);
elseif(bitCoding == 16)
    finalResult = uint16(finalResult);
else
    disp('Image is either 8 or 16-bit, you fool!');
end

%save the result
imwrite(finalResult, 'results/LF/_init_nocorr.png');


%parallel part processing all files
parfor k = 1 : length(lfield)
    %needs to be re-set for some reason otherwise Matlab does not compile this file
    colourSpace = varargin{3};
    
    fileName = lfield(k).name;
    he1 = imread(fileName);
    
    % if(bitCoding == 8)
    %     he1 = im2uint8(he1);
    % end
    
    if(strcmp(colourSpace, 'CIELab') && strcmp(clusterFun, 'KMeans'))
        he1	= rgb2lab(he1);
    elseif(strcmp(colourSpace, 'CIELab') && strcmp(clusterFun, 'MVQ'))
        disp('MVQ clustering cannot be applied in CIELAB space in this implementation. Using RGB space instead.');
        colourSpace = 'RGB';
    end

    dfull2 = double(he1);
    fnrows2 = size(dfull2,1);
    fncols2 = size(dfull2,2);
    full_transform2 = reshape(dfull2,fnrows2*fncols2,3);
    
    fprintf(1, 'Applying colour transfer to image : %s\n', fileName);
    fullTransform2 = mg_transform_tps_parallel(param, full_transform2, config.ctrl_pts);
    disp('Finished.');

    if(strcmp(colourSpace, 'CIELab'))
        fullTransform2 = lab2rgb(fullTransform2);
        ind3 = find(fullTransform2 < 0);
        fullTransform2(ind3) = 0;
        ind4 = find(fullTransform2 > 1);
        fullTransform2(ind4) = 1;
        fullTransform2 = 255*fullTransform2;
    end
    
    finalResult2 = reshape(fullTransform2, [fnrows2 fncols2 3]);
    if(bitCoding == 8)
        finalResult2 = uint8(finalResult2);
    elseif(bitCoding == 16)
        finalResult2 = uint16(finalResult2);
    else
        disp('Image is either 8 or 16-bit, you fool!');
    end

    %save the result
    outputFileName = strcat('results/LF/nocorr_', fileName);
    imwrite(finalResult2, outputFileName);
end

toc
end

% call example :
% ctfunction_parallel('data/strelitzia.jpg', 'data/eucalyptus/', 'data/eucalyptus/in_08_08_3012.329346_1538.333984_.png', 8, 'MVQ', 50, 'RGB')