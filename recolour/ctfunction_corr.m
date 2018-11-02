function [param, config] = ctfunction_corr(palette, target, bitCoding, annsteps, colPalette, colTarget, outname, param_init, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ctfunction_corr(paletteFile, targetFile, colourSpace, numCorr)
% Computes the colour transfer result between a target image and palette
% image that are aligned, eg. pixels at the same location in both images
% are correspondences.
%
% palette:      Path to palette file to be processed. The colour
%               distribution of the palette file will be mapped to the target image. 
% target:       Path to target file to be processed.
% bitCoding:    Allows to process 16-bit images instead of the 8-bit regularly.
% annstep:      Change the amount of steps in the annealing process used for the optimisation.
% colPalette:   Colours corresponding to the matched pixels in im1(target)
% colTarget:    Colours corresponding to the matched pixels in im2(palette)
% outname:      File name for the output image
% param_init:   Parameters for the colour transfer (default or precomputed)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin < 8)
    % Eight input arguments must be given
    error('Eight input arguments are required: palette, target, bit encoding, anneal steps, matchcol of pal, matchcol of tar, result name, initial param');
elseif(nargin == 8)
    % Optionnally you can specify the colour space you want to work in
    colourSpace = 'RGB';
elseif(nargin == 9)
    colourSpace = varargin{1};
end
%addpath(genpath('../L2RegistrationForCT'));

close all;

%read target image and reshape
he1 = target;
if(strcmp(colourSpace, 'CIELab'))
    colTarget = colTarget./255;
    colTarget = rgb2lab(colTarget);
    colPalette = colPalette./255;
    colPalette = rgb2lab(colPalette);
    he1	= rgb2lab(he1);
end
dfull1 = double(he1);
fnrows = size(dfull1,1);
fncols = size(dfull1,2);
ab1 = reshape(dfull1,fnrows*fncols,3);


he2 = palette;
if(bitCoding == 16)
    he2 = im2uint16(he2);
end
if(strcmp(colourSpace, 'CIELab'))
	he2	= rgb2lab(he2);
end
%ab2 = double(he2);
%nrows = size(ab2,1);
%ncols = size(ab2,2);
%ab2 = reshape(ab2,nrows*ncols,3);
%ind = randperm(nrows*ncols,3*round((nrows*ncols)/100));
%ind = randperm(nrows*ncols,numCorr);
%ind = randperm(min(fnrows,nrows)*min(fncols,ncols),numCorr);
%X = colX; Y = colY;

[idx,colMean] = kmeans(colTarget,1000); %apply kmeans to find 1000 even sampled colours to get good initial guess of parameters
IdTar = knnsearch(colTarget,colMean); %find the correspondences closest to k means
corrTarget = colTarget(IdTar,:);
corrPalette = colPalette(IdTar,:);
[n,d] = size(corrPalette);

%initialise the parameters for registration 
[config] = mg_initialize_config_corr(corrTarget, corrPalette, colourSpace, bitCoding, annsteps, param_init);
disp('Registration of colours started...');

%begin with rough estimate of parameters using 1000 colours
config.max_iter = 500; %limit iterations to 500 
for i = 1:config.AnnSteps
    config.iter = (config.AnnSteps-i+1);
    if(i == config.AnnSteps )
        config.max_iter = 800; %increase iterations at final annealling step
    end
    [param, transformed_model, history, config] = gmmreg_L2_corr(config, 0:(n-1), 0:(n-1));
    if(i ~= config.AnnSteps)
        config.scale = .5*config.scale;
    end
end

%estimate more accurate parameters using all colours correspondences
config.model = colTarget;
config.scene = colPalette;
[n,d] = size(colPalette);
config.max_iter = 10000; %set max iter 
config.opt_affine = 0; %we do not estimate the affine part in this step
[param, transformed_model, history, config] = gmmreg_L2_corr(config, 0:(n-1), 0:(n-1));

disp('Registration of colours finished.');

%apply the colour transformation to the target image to get the result
%image
disp('Applying colour transfer to target image...');
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
%imshow(finalResult);

%save the result
%disp("Saving the result in results/corrresult.png"); 

imwrite(finalResult, outname);

end