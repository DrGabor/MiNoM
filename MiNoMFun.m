%%%%%%%% MiNoM scan matching algorithm. 
function [rstR, rstT, TModel] = MiNoMFun(params)
if nargin == 0
    clc; close all; clear all;
    DataDir = 'obs_000000.ply'; 
    cloud_mov = pcread(DataDir);
    DataDir = 'obs_000002.ply';
    cloud_ref = pcread(DataDir);
    params = genParamsFun(cloud_mov.Location', cloud_ref.Location', ...
            'ref_normal', cloud_ref.Normal', ...
            'P', [1.0 2.0], ...
            'mode', 'point2plane', ...
            'optimize', 'IRLS', ...
            'is_show', 1, ...
            'verbose', 1);
end
Mov0 = params.Mov;
Ref0 = params.Ref;
params.Ref0 = Ref0; 

nDim = size(Mov0, 1);
P = params.P;
K = length(P);
R0 = params.Tf0(1:nDim, 1:nDim);
t0 = params.Tf0(1:nDim, end);
%%%%% transform the Mov point cloud. 
Mov0 = Loc2Glo(Mov0, R0', t0);
Mov_Normal = params.Mov_Normal; 
if ~isempty(Mov_Normal)
    Mov_Normal = R0 * Mov_Normal;
end
params.Mov0 = Mov0; 
params.Mov_Normal = Mov_Normal; 

LLH = [];
MdRef = createns(Ref0');
nCnt = 0;
TModel = [];
tmp = [];
tmp.tf = eye(nDim+1); % params.Tf0(1:nDim, :);
TfArray = tmp;
if strcmpi(params.mode, 'point2plane')
    Ref_Normal0 = params.Ref_Normal;
end
if strcmpi(params.mode, 'plane2plane')
    kNum = 20;
    Cov_Mov0 = compute_covariance(Mov0, kNum);
    Cov_Ref0 = compute_covariance(Ref0, kNum);
end
max_irls_iter = 1; 
for nIter_out = 1 : 1 : params.maxIter_icp
    dTf = TfArray(end).tf;
    dR = dTf(1:nDim, 1:nDim);
    dT = dTf(1:nDim, end);
    %% corespondence, ratio verification.
    Aft = Loc2Glo(Mov0, dR', dT);
    [NNIdx, DD] = knnsearch(MdRef, Aft');
    EffIdx = find(DD <= inf);
    Ref = Ref0(:, NNIdx(EffIdx));
    Aft = Aft(:, EffIdx);
    Mov = Mov0(:, EffIdx);
    params.Mov = Mov;
    params.Ref = Ref;
    params.Aft = Aft;
    params.Tf  = dTf;
    if strcmpi(params.mode, 'point2plane')
        params.Ref_Normal = Ref_Normal0(:, NNIdx(EffIdx));
    end
    [Err, MArray] = compute_residualFun(params);
    params.MArray = MArray;
    %%
    [Model, Z] = MoEPFittingMultiFun(Err, P, params.maxEMCounter, params.maxIter_em, params.bRand, 0);
    LLH(end+1) = llhFun( Err, Model );
    if strcmpi(params.optimize, 'irls')
        
        [dR, dT] = IRLSFun( params, Z, Model, max_irls_iter );
        max_irls_iter = min( max_irls_iter + 1, params.maxIter_irls ); 
    end
    R = TfArray(end).tf(1:nDim, 1:nDim);
    T = TfArray(end).tf(1:nDim, end);
    R = dR * R;
    T = dR * T + dT;
    tmp = [];
    tmp.tf = [R T];
    TfArray(:, end+1) = tmp;
    %% terminating condition.
    if length(LLH) >= 2
        if abs( LLH(end) - LLH(end-1) ) / abs(LLH(end-1) ) <= params.tf_eps % 1e-3
            nCnt = nCnt + 1;
        end
    end
    if nCnt >= 3
        break;
    end
    tmp = [];
    tmp.model = Model;
    tmp.Tf = [R T];
    tmp.llh = LLH(end);
    TModel = [TModel tmp];
    str = sprintf('nIter = %03d, mode = %s, llh = %.6f', nIter_out, params.mode, LLH(end) );
    if params.verbose
        disp(str); 
        % disp([R T]);
    end
end
if params.IS_SHOW
    VisualizeFun(TModel, params);
end
opt_mix = TfArray(end).tf;
dR = opt_mix(1:nDim, 1:nDim);
dT = opt_mix(1:nDim, end);
Tf = [dR dT; zeros(1, nDim) 1] * [R0 t0; zeros(1, nDim) 1];
rstR = Tf(1:nDim, 1:nDim);
rstT = Tf(1:nDim, end);
end
