function [config] = mg_initialize_config_corr(model, scene, colourSpace, bitCoding, annsteps, varargin)
config.model = model;
config.scene = scene;
% estimate the scale from the covariance matrix
[n,d] = size(model);
config.scale = 13;%power(det(model'*model/n), 1/(2^d));
config.display = 0;
config.init_param = [ ];
config.max_iter = 25000;
config.normalize = 0;
config.functionType = 'TPS';
config.AnnSteps = annsteps;
config.scale = (2^(config.AnnSteps-1))*(config.scale);


switch lower(config.functionType)
    case 'tps'
        interval = 5;%found that the best value is 5. 
        if(bitCoding == 8)
            config.ctrl_pts =  set_ctrl_pts8(model, scene, interval, d, colourSpace);%set control points in a regular grid spanning the colour space 
        elseif(bitCoding == 16)
            config.ctrl_pts =  set_ctrl_pts16(model, scene, interval, d, colourSpace);
        end
        config.beta = 0.003; % this value controls the strength of the regularisation term, we found 0.003 gives the best results when correspondences used
        config.alpha = 1 - config.beta;
        config.opt_affine = 1;
        [n,d] = size(config.ctrl_pts); % number of points in model set
        if(nargin == 6 && ~isempty(varargin{1}))
            param_init = varargin{1};
            config.init_param = param_init;
            config.init_tps = param_init(end+1-d*(n-d-1):end);
            config.init_affine = param_init(1:d*(d+1));
        else
            config.init_tps = zeros(n-d-1,d);
            init_affine = repmat([zeros(1,d) 1],1,d);
            config.init_param = [init_affine zeros(1, d*n-d*(d+1))];
            config.init_affine = [ ];
        end
end


