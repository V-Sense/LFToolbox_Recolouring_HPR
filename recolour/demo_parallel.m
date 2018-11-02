function [] = demo_parallel(paletteFile, newPaletteFile, targetFolder, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo_parallel(paletteFile, newPaletteFile, targetFolder, varargin)
% Computes the colour transfer result between a target image and palette
% image, then applies the colour transfer using the first target (now
% recoloured) to every image in a light field array
%
% paletteFile:      The name of the palette file to be processed. The colour
%                   distribution of the palette file will be mapped to the target image. 
% newPaletteFile:   The name of the target file to be processed. Will be later
%                   used in the recolouring of the rest of the LF array.
% targetFolder:     The name of the folder containing the LF array images.
% clusterFun:       'KMeans' or 'MVQ'(default). The clustering function used.
% nColors:          The number of clusters computed by the clustering
%                   function. The default is 50.
% colourSpace:      The colour space in which the registration is performed. The default is 'RGB'. 
% numCorr:          The number of correspondences used to compute the colour transfer function.
%                   Default is 5000.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(nargin < 3)
    %Three input arguements must be given
    error('Two input arguements are required: palette image, target image.');
elseif(nargin == 3)
    %If a clustering function is not specified, use Matlab's MVQ fucntion.
    %If the number of clusters to be used is not specified, use 50.
    %If the colour space is not specified, use RGB.
    %If the number of correspondences to be used is not specified, use 5000.
    clusterFun = 'MVQ';
    nColors = 50;
    colourSpace = 'RGB';
    numCorr = 5000;
elseif(nargin == 4)
    clusterFun = varargin{1};
    nColors = 50;
    colourSpace = 'RGB';
    numCorr = 5000;
elseif(nargin == 5)
    clusterFun = varargin{1};
    nColors = varargin{2};
    colourSpace = 'RGB';
    numCorr = 5000;
elseif(nargin == 6)
    clusterFun = varargin{1};
    nColors = varargin{2};
    colourSpace = varargin{3};
    numCorr = 5000;
elseif(nargin == 7)
    clusterFun = varargin{1};
    nColors = varargin{2};
    colourSpace = varargin{3};
    numCorr = varargin{4};
end

test = imread(newPaletteFile);
disp("showing image maybe perhaps..");
imshow(test);

ctfunction(paletteFile, newPaletteFile, clusterFun, nColors, colourSpace);

ctfunction_corr_parallel('results/LF/_recolour.png', targetFolder, colourSpace, numCorr);

%call example : demo_parallel('data/isaac.png', 'data/eucalyptus/out_08_08_3012.329346_1538.333984_.png', 'data/eucalyptus', 'MVQ', 70, 'RGB', 5000)
%For 'best' results, try to pick a newPaletteFile (arg 2) that is situated roughly in the middle of the LF array.