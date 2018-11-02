function [param, config] = ctfunction_corr_parallel_all(paletteFile, targetFolder, bitCoding, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ctfunction_corr(paletteFile, targetFile, colourSpace, numCorr)
% Computes the colour transfer result between a target image and palette
% image that are aligned, eg. pixels at the same location in both images
% are correspondences.
%
% paletteFile:  The name of the palette file to be processed. The colour
%               distribution of the palette file will be mapped to the target image. 
% targetFile:   The name of the target file to be processed.The colour
%               distribution of the palette file will be mapped to the target image.
% colourSpace:  'RGB' or 'CIELab' . The colour space in which the registration is performed. 
% numCorr:      The number of correspondences used to compute the colour transfer function, eg. 50
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if(nargin < 3)
    %Three input arguements must be given
    error('Three input arguements are required: palette image, target folder, target image bit encoding');
elseif(nargin == 3)
    colourSpace = 'RGB';
    numCorr = 5000;
elseif(nargin == 4)
    colourSpace = varargin{1};
    numCorr = 5000;
elseif(nargin == 5)
    colourSpace = varargin{1};
    numCorr = varargin{2};
end
addpath(genpath('../L2RegistrationForCT'));

close all;

LFfolder = targetFolder;

if ~isdir (LFfolder);
    errorMessage = sprintf('Error : following folder does not exist in this realm : \n%s', LFfolder);
    uiwait(warndlg(errorMessage));
    return;
end

filePattern = fullfile(LFfolder, '*.png');
lfield = dir(filePattern);


he2 = imread(paletteFile);
if(bitCoding == 16)
    he2 = im2uint16(he2);
end
if(strcmp(colourSpace, 'CIELab'))
	he2	= rgb2lab(he2);
end
ab2 = double(he2);
nrows = size(ab2,1);
ncols = size(ab2,2);
ab2 = reshape(ab2,nrows*ncols,3);

%imshow(finalResult);
% 
% %save the result
% %disp("Saving the result in Results/corrresult.png"); 
% 
% outputFileName = strcat('Results/LF/', '_recolour.png');
% imwrite(finalResult, outputFileName);

parfor k = 1 : length(lfield)
    fileName = lfield(k).name;
    he1 = imread(fileName);
    if(bitCoding == 8)
        he1 = im2uint8(he1);
    end
    if(strcmp(colourSpace, 'CIELab'))
        he1 = rgb2lab(he1);
    end

    dfull1 = double(he1);
    fnrows = size(dfull1,1);
    fncols = size(dfull1,2);
    ab1 = reshape(dfull1,fnrows*fncols,3);
    
    ind = randperm(min(fnrows,nrows)*min(fncols,ncols),numCorr);
    X = ab1(ind,:); Y = ab2(ind,:);
    [n,d] = size(X);

    %initialise some parameters used to control the registration 
    [config] = mg_initialize_config_corr(X,Y, colourSpace, bitCoding);
    disp('Registration of colours started...');
    %register the colours using an annealling scheme, estimate a TPS transformation 
    for i = 1:config.AnnSteps
        config.iter = (config.AnnSteps-i+1);
        [param, transformed_model, history, config] = gmmreg_L2_corr(config, 0:(n-1), 0:(n-1));
        if(i ~= config.AnnSteps)
            config.scale = .5*config.scale;
        end
    end
    disp('Registration of colours finished.');
    
    fprintf(1, 'Applying colour transfer to image : %s\n', fileName);
    fullTransform = mg_transform_tps_parallel(param, ab1, config.ctrl_pts);
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
    
    outputFileName = strcat('results/LF/reco_corr_', fileName);
    imwrite(finalResult, outputFileName);
end

toc
end

%call example : ctfunction_corr_parallel_all('data/artichoke.jpg', 'data/eucalyptus/', 8, 'RGB', 5000)