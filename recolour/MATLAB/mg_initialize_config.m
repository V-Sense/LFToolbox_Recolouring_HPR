function [config] = mg_initialize_config(model, scene, colourSpace, bitCoding)

config.model = model;
config.scene = scene;
% estimate the scale from the covariance matrix
[n,d] = size(model);
config.scale = power(det(model'*model/n), 1/(2^d));
config.display = 0;
config.init_param = [ ];
config.max_iter = 10000;
config.normalize = 0;
config.functionType = 'tps';
if(bitCoding == 8)
    config.AnnSteps = 5; %found that the best value it 5. 
elseif(bitCoding == 16)
    config.AnnSteps = 8; %who knows what the best value is for 16bit images..
end
config.scale = (2^(config.AnnSteps-1))*(config.scale);
switch lower(config.functionType)
    case 'tps'
        interval = 5;%found that the best value it 5.
        if(bitCoding == 8)
            config.ctrl_pts =  set_ctrl_pts8(model, scene, interval, d, colourSpace);%set control points in a regular grid spanning the colour space 
        elseif(bitCoding == 16)
            config.ctrl_pts =  set_ctrl_pts16(model, scene, interval, d, colourSpace);
        end            
        config.beta = 0.000003; % this value controls the strength of the regularisation term, we found 0.000003 gives the best results. 
        config.alpha = 1 - config.beta;
        config.opt_affine = 1;
        [n,d] = size(config.ctrl_pts); % number of points in model set
        config.init_tps = zeros(n-d-1,d);
        init_affine = repmat([zeros(1,d) 1],1,d);
        config.init_param = [init_affine zeros(1, d*n-d*(d+1))];
        config.init_affine = [ ];
end

