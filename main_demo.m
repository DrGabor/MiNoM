clc; close all;
DataDir = 'raw_000802.pcd';
cloud_mov = pcread(DataDir);
cloud_mov = pcdownsample(cloud_mov, 'random', 0.1);  % downsample for efficiency. 
DataDir = 'raw_000804.pcd';
cloud_ref = pcread(DataDir);
params = genParamsFun(cloud_mov.Location', ...  % 3XN matrix.
     cloud_ref.Location', ...                   % 3XN matrix.
    'ref_normal', cloud_ref.Normal', ...        % 3XN matrix. 
    'P', [1.0 2.0], ...        % [2.0 2.0] can also be tried, but [1.0 2.0] seems to be more robust. 
    'mode', 'point2plane', ... % 'point2plane' is robust than 'point2point' 
    'Tf0', eye(4), ...
    'is_show', 1, ... % is_show = 0 when use MiNoM in other applications.  
    'verbose', 1 );   % verbose = 0 when use MiNoM in other applications. 
[dR, dT] = MiNoMFun(params);
