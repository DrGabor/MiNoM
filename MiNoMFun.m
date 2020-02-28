%%%%%%%% MiNoM scan matching algorithm. 
function [rstR, rstT, TModel] = MiNoMFun(params)
if nargin == 0
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
params.Mov0 = Mov0; 

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
    [Model, Z] = MoEPFittingFun(Err, P, params.maxIter_em, params.bRand);
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

function [varargout] = compute_residualFun(params)
Err = []; 
MArray = []; 
Aft = params.Aft;
Ref = params.Ref;
V = Aft - Ref; 
if strcmpi(params.mode, 'point2point')
    Err = sqrt(sum(V.^2));
end
if strcmpi(params.mode, 'point2plane')
    n = params.Ref_Normal; 
    tmp = V .* n;
    Err = abs(sum(tmp)); % DD(EffIdx)';
    if params.Hb_point2plane   
        MArray = [n(1, :).^2
            n(1, :).*n(2, :)
            n(1, :).*n(3, :)
            n(2, :).^2
            n(2, :).*n(3, :)
            n(3, :).^2];
    end
end
if nargout == 1
    varargout{1} = Err; 
end
if nargout == 2
    varargout{1} = Err; 
    varargout{2} = MArray; 
end
end

function [varargout] = MoEPFittingFun(Err, P0, maxIter_em, bRand)
%%%%%% make sure that each time produces same random number.
% rng('default');
% rng(10000);
if isrow(P0)
    P0 = P0';
end
if isrow(Err)
    Err = Err';
end
ss = struct('type', [], 'para', [], 'prior', [], 'coeff', []);
Model = repmat(ss, 1, length(P0));
K = length(Model);
for k = 1 : 1 : K
    Model(k).type = 'ep';
    Model(k).para = 1.0; % rand;
    Model(k).prior = P0(k);
end
if bRand
    Coeff = rand(K, 1); 
else
    Coeff = ones(K, 1);  %   %
end
Coeff = Coeff / sum(Coeff);
for id = 1 : 1 : K
    Model(id).coeff = Coeff(id);
end
% llh_Loc = [];
for i = 1 : 1 : maxIter_em
    %%%%% given a error, estimate paramters.
    % E step
    Z = expectationFun( Err, Model );
    % M step
    Model = maximationFun( Err, Z, Model );
    % llh_Loc(end+1) = llhFun( Err, Model );
end
varargout{1} = Model;
if nargout == 2
    varargout{2} = Z;
end
end

function Z = expectationFun( Err, Model )
K = length(Model);
N = length(Err);
Prob = zeros(K, N);
for k = 1 : 1 : K
    if strcmp(Model(k).type, 'ep')
        p = Model(k).prior;
        sita = Model(k).para;
        pp = p*sita^(1/p) / (2*gamma(1/p)) * exp(-sita*abs(Err).^p);
    end
    if strcmp(Model(k).type, 'uniform')
        pp = Model(k).para * ones(1, N);
    end
    Prob(k, :) = Model(k).coeff * pp;
end
Prob_Sum = sum(Prob, 1);
Prob_Sum(Prob_Sum < 1e-8) = 1e-8;
Z = zeros(K, N);
for k = 1 : 1 : K
    Z(k, :) = Prob(k, :) ./ Prob_Sum;
end
end

function Model = maximationFun( Err, Z, Model0 )
if ~isrow(Err)
    Err = Err';
end
Model = Model0;
for nIter = 1 : 1 : 1
    %% estimate parameters of EP.
    [Model, Omega] = weightFun(Z, Model);
    for k = 1 : 1 : length(Model)
        if strcmp(Model(k).type, 'ep')
            p = Model(k).prior;
            tmp = Z(k, :) .* Err.^p;
            Model(k).para = Omega(k) / (p * sum(tmp));
        end
    end
end
end

function [ModelNew, Omega] = weightFun(Z, Model)
Omega = sum(Z, 2);
tmp = Omega/sum(Omega); % ( - Lamda*Df)/(1-Lamda*K*Df);
Coeff = max(0, tmp);
Coeff = Coeff / sum(Coeff);
ModelNew = Model;
for k = 1 : 1 : length(Model)
    ModelNew(k).coeff = Coeff(k);
end
end

function [ llh ] = llhFun( data, Model )
if ~isrow(data)
    data = data';
end
Dim = size(data, 1);
N = size(data, 2);
K = length(Model);

Prob = zeros(K, N);
%%%%%% estimate posterior.
for k = 1 : 1 : K
    if strcmp(Model(k).type, 'ep')
        sita = Model(k).para;
        p = Model(k).prior;
        pp = 2*eppdfFun(data, sita, p);
    end
    if strcmp( Model(k).type, 'uniform')
        Prob(k, :) = Model(k).para * ones(1, N);
    end
    Prob(k, :) = Model(k).coeff * pp;
end
Prob_Sum = sum(Prob);
Prob_Sum(Prob_Sum < 1e-8) = 1e-8;
llh = sum(log(Prob_Sum));
llh = llh / length(Prob_Sum);
end

function [ y ] = eppdfFun( data, sita, p )
lamda = p*sita^(1/p)/(2*gamma(1/p));
y = lamda * exp(-sita*(abs(data)).^p);
end

function [ varargout ] = IRLSFun( params, Z0, Model, maxIter )
Z = [];
Sita = [];
P = [];
for k = 1 : 1 : length(Model)
    if strcmp(Model(k).type, 'ep')
        Z = [Z; Z0(k, :)];
        Sita = [Sita; Model(k).para];
        P = [P Model(k).prior];
    end
end
N = size(Z, 2);
Dim = size(params.Ref, 1);
% W = Z .* repmat(Sita, 1, N);
W = ones( 1, length(params.Ref) );
%%%%%%%%% optimize Aft + dX - Ref alternatively.
W = ones(1, N);
params.W = W;
Aft = params.Aft;
Ref = params.Ref;
XArray = [];
for i = 1 : 1 : maxIter
    %% update transformation.
    if strcmpi(params.mode, 'point2point')
        [dR, dT] = reg_point2point(params);
    end
    if strcmpi(params.mode, 'point2plane')
        [dR, dT] = reg_point2plane(params);
    end
    if strcmpi(params.mode, 'plane2plane')
        [dR, dT] = reg_plane2plane(params);
    end
    tmp = [];
    tmp.Tf = [dR dT];
    XArray = [XArray tmp];
    %% update W.
    tmpAft = Loc2Glo(Aft, dR', dT);
    % tmpAft = bsxfun(@plus, dX, Aft);
    tmpDiff = tmpAft - Ref;
    r = sqrt( sum(tmpDiff.^2));
    r( r < 1e-8) = 1e-8;  % regularition for numerical stable.
    W = zeros(1, N);
    for k = 1 : 1 : size(Z, 1)
        p = P(k);
        W = W + Z(k, :) * Sita(k) .* r .^(p-2);
    end
    W = W/sum(W);
    params.W = W;
end
if nargout == 1
    varargout{1} =  XArray;
end
if nargout == 2
    varargout{1} = dR;
    varargout{2} = dT;
end
end

function [dR, dT] = reg_point2plane(params)
if params.Hb_point2plane
    [H, a] = CalHbFun(params, 1);
    dx = -pinv(H)*a;
    dR = expm(SkewFun(dx(1:3)));
    Tf = params.Tf;
    R0 = Tf(1:3, 1:3);
    dT = R0'*dx(4:6);
else
    s = params.Aft;
    d = params.Ref;
    n = params.Ref_Normal;
    W = params.W;
    V = d - s;
    a = sum( n .* V)';
    Dim = size(s, 1);
    
    if Dim == 3
        a1 = n(3, :) .* s(2, :) - n(2, :) .* s(3, :);
        a2 = n(1, :) .* s(3, :) - n(3, :) .* s(1, :);
        a3 = n(2, :) .* s(1, :) - n(1, :) .* s(2, :);
        A = [a1' a2' a3' n'];
    else
        a1 = n(2, :).*s(1, :) - n(1, :).*s(2, :);
        A = [a1' n'];
    end
    if ~isrow(W)
        W = W';
    end
    W2 = sqrt(W);
    A = repmat(W2', 1, size(A, 2)) .* A;
    a = W2' .* a;
    %%%%%%%%%%% solve the problem || Ax - b ||, the solution also can be x = pinv( A' * A) * A' * b - x
    H = A'*A;
    b = -A'*a;
    x = -pinv(H) * b;
    % [dR, dT] = ObtainTf(x);
    Ang = x(3:-1:1);
    dR = eul2rotm(Ang');
    dT = x(4:6);
end
end

%%%%%%%%%%%%%%% Calculate the registration matrix %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% T(TData)->MData %%%%%%%%%%%%%%%%%%%%%%%%%
% SVD solution
function [R1, t1] = reg_point2point(varargin)
if nargin == 1
    params = varargin{1}; 
end
if nargin == 3
    params = []; 
    params.Aft = varargin{1}; 
    params.Ref = varargin{2}; 
    params.W   = varargin{3}; 
end
MovData = params.Aft;
RefData = params.Ref;
W       = params.W;
%%%%%%% normalize the weight.
Dim = size(MovData, 1);
W = W / sum(W);
W_Normalize = repmat(W, Dim, 1);
M = RefData;
mm = sum(M.*W_Normalize, 2);
S = MovData;
ms = sum(S.*W_Normalize, 2);
Sshifted = bsxfun(@minus, S, ms );
Mshifted = bsxfun(@minus, M, mm );

Sshifted = Sshifted .* W_Normalize;
K = Sshifted*Mshifted';
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
    B = eye(Dim);
    B(Dim,Dim) = det(V*U');
    R1 = V*B*U';
end
t1 = mm - R1*ms;
end

function [H, b] = CalHbFun(params, bFast)
Aft = params.Aft;
Ref = params.Ref;
Mov = params.Mov;
w   = params.W;
w = w / sum(w);
MArray = params.MArray;
Tf = params.Tf;
R = Tf(1:3, 1:3);
PtsDiff = Aft - Ref;
if bFast
    %%%%%% assign to variables.
    m1 = MArray(1,:);
    m2 = MArray(2,:);
    m3 = MArray(3,:);
    m4 = MArray(4,:);
    m5 = MArray(5,:);
    m6 = MArray(6,:);
    R1_1 = R(1);
    R1_2 = R(4);
    R1_3 = R(7);
    R2_1 = R(2);
    R2_2 = R(5);
    R2_3 = R(8);
    R3_1 = R(3);
    R3_2 = R(6);
    R3_3 = R(9);
    v1 = PtsDiff(1,:);
    v2 = PtsDiff(2,:);
    v3 = PtsDiff(3,:);
    p1 = Mov(1,:);
    p2 = Mov(2,:);
    p3 = Mov(3,:);
    t2 = R1_2.*p3;
    t6 = R1_3.*p2;
    t3 = t2-t6;
    t4 = R2_2.*p3;
    t7 = R2_3.*p2;
    t5 = t4-t7;
    t8 = R3_2.*p3;
    t10 = R3_3.*p2;
    t9 = t8-t10;
    t11 = m1.*t3.*w;
    t12 = m2.*t5.*w;
    t13 = m3.*t9.*w;
    t14 = t11+t12+t13;
    t15 = m2.*t3.*w;
    t16 = m4.*t5.*w;
    t17 = m5.*t9.*w;
    t18 = t15+t16+t17;
    t19 = m3.*t3.*w;
    t20 = m5.*t5.*w;
    t21 = m6.*t9.*w;
    t22 = t19+t20+t21;
    t23 = R1_1.*p3;
    t29 = R1_3.*p1;
    t24 = t23-t29;
    t25 = R2_1.*p3;
    t30 = R2_3.*p1;
    t26 = t25-t30;
    t27 = R3_1.*p3;
    t31 = R3_3.*p1;
    t28 = t27-t31;
    t32 = m1.*t24.*w;
    t33 = m2.*t26.*w;
    t34 = m3.*t28.*w;
    t35 = t32+t33+t34;
    t36 = m2.*t24.*w;
    t37 = m4.*t26.*w;
    t38 = m5.*t28.*w;
    t39 = t36+t37+t38;
    t40 = m3.*t24.*w;
    t41 = m5.*t26.*w;
    t42 = m6.*t28.*w;
    t43 = t40+t41+t42;
    t44 = R1_1.*p2;
    t50 = R1_2.*p1;
    t45 = t44-t50;
    t46 = R2_1.*p2;
    t51 = R2_2.*p1;
    t47 = t46-t51;
    t48 = R3_1.*p2;
    t52 = R3_2.*p1;
    t49 = t48-t52;
    t53 = m1.*t45.*w;
    t54 = m2.*t47.*w;
    t55 = m3.*t49.*w;
    t56 = t53+t54+t55;
    t57 = m2.*t45.*w;
    t58 = m4.*t47.*w;
    t59 = m5.*t49.*w;
    t60 = t57+t58+t59;
    t61 = m3.*t45.*w;
    t62 = m5.*t47.*w;
    t63 = m6.*t49.*w;
    t64 = t61+t62+t63;
    t65 = -t11-t12-t13;
    t66 = -t53-t54-t55;
    t67 = -t15-t16-t17;
    t68 = -t57-t58-t59;
    t69 = m2.*w;
    t70 = -t19-t20-t21;
    t71 = -t61-t62-t63;
    t72 = m3.*w;
    t73 = m5.*w;
    
    A0(0+1,0+1) = sum(t3.*t14+t5.*t18+t9.*t22);
    A0(0+1,1+1) = sum(-t14.*t24-t18.*t26-t22.*t28);
    A0(0+1,2+1) =  sum(t14.*t45+t18.*t47+t22.*t49);
    A0(0+1,3+1) =  sum(t65);
    A0(0+1,4+1) =  sum(t67);
    A0(0+1,5+1) =  sum(t70);
    A0(0+1,6+1) =  sum(-t14.*v1-t18.*v2-t22.*v3);
    A0(1+1,0+1) =  sum(-t3.*t35-t5.*t39-t9.*t43);
    A0(1+1,1+1) =  sum(t24.*t35+t26.*t39+t28.*t43);
    A0(1+1,2+1) =  sum(-t35.*t45-t39.*t47-t43.*t49);
    A0(1+1,3+1) =  sum(t35);
    A0(1+1,4+1) =  sum(t39);
    A0(1+1,5+1) =  sum(t43);
    A0(1+1,6+1) =  sum(t35.*v1+t39.*v2+t43.*v3);
    A0(2+1,0+1) =  sum(t3.*t56+t5.*t60+t9.*t64);
    A0(2+1,1+1) =  sum(-t24.*t56-t26.*t60-t28.*t64);
    A0(2+1,2+1) =  sum(t45.*t56+t47.*t60+t49.*t64);
    A0(2+1,3+1) =  sum(t66);
    A0(2+1,4+1) =  sum(t68);
    A0(2+1,5+1) =  sum(t71);
    A0(2+1,6+1) =  sum(-t56.*v1-t60.*v2-t64.*v3);
    A0(3+1,0+1) =  sum(t65);
    A0(3+1,1+1) =  sum(t35);
    A0(3+1,2+1) =  sum(t66);
    A0(3+1,3+1) =  sum(m1.*w);
    A0(3+1,4+1) =  sum(t69);
    A0(3+1,5+1) =  sum(t72);
    A0(3+1,6+1) =  sum(w.*(m1.*v1+m2.*v2+m3.*v3));
    A0(4+1,0+1) =  sum(t67);
    A0(4+1,1+1) =  sum(t39);
    A0(4+1,2+1) =  sum(t68);
    A0(4+1,3+1) =  sum(t69);
    A0(4+1,4+1) =  sum(m4.*w);
    A0(4+1,5+1) =  sum(t73);
    A0(4+1,6+1) =  sum(w.*(m2.*v1+m4.*v2+m5.*v3));
    A0(5+1,0+1) =  sum(t70);
    A0(5+1,1+1) =  sum(t43);
    A0(5+1,2+1) =  sum(t71);
    A0(5+1,3+1) =  sum(t72);
    A0(5+1,4+1) =  sum(t73);
    A0(5+1,5+1) =  sum(m6.*w);
    A0(5+1,6+1) =  sum(w.*(m3.*v1+m5.*v2+m6.*v3));
    H = A0(:, 1:6);
    b = A0(:, end);
else
    H = zeros(6, 6);
    b = zeros(6, 1);
    parfor id = 1 : 1 : size(PtsDiff, 2)
        v = PtsDiff(:, id);
        pt_mov = Mov(:, id);
        A = [-R*SkewFun(pt_mov) eye(3)];
        w = W(id);
        m = MArray(:, id);
        M = [m(1) m(2) m(3)
            m(2) m(4) m(5)
            m(3) m(5) m(6)];
        H = H + w * A'*M*A;
        b = b + w * A'*M*v;
    end
end
end

function VisualizeFun(TModel, params)
Mov0 = params.Mov0; 
Ref0 = params.Ref0; 
Aft = params.Aft;
optModel = TModel(end).model;
tmp = struct('data', []);
P = cat(1, optModel(:).prior);
TSita = repmat(tmp, 1, length(P));
for id = 1 : 1 : length(TModel)
    optModel = TModel(id).model;
    for k = 1 : 1 : length(optModel)
        if strcmp(optModel(k).type, 'ep')
            p = optModel(k).prior;
            [~, ind] = ismember(p, P);
            TSita(ind).data(end+1) = optModel(k).para;
        end
    end
end
h = figure;
subplot(121);
hold on;
grid on;
llh = cat(1, TModel(:).llh);
plot(llh, 'b.-');
title('Likelihood Function');
xlabel('Iteration');
ylabel('llh');
%%%%%%%%%%
figure(h);
subplot(122);
hold on;
grid on;
Str = {};
h = [];
ColorArray = {'r', 'b', 'm', 'k', 'g', 'c'};
for id = 1 : 1 : length(TSita)
    y = TSita(id).data;
    colorId = mod(id, length(ColorArray));
    if colorId == 0
        colorId = length(ColorArray);
    end
    color = ColorArray{colorId};
    h = [h plot(y, 'marker', '.', 'linestyle', '-', 'color', color)];
    Str{end+1} = sprintf('p = %.1f', P(id));
end
xlabel('Iteration');
ylabel('$\bf{\theta}$ value', 'interpreter', 'latex');
legend(h, Str);
legend('boxoff');
legend('location', 'best');
str = 'covergence curve for $\bf{\theta}$';
title(str, 'interpreter', 'latex');
%%%%%%%%%%
%% plot registration results.
nDim = size(Ref0, 1);
figure;
hold on;
grid on;
% axis equal;
if nDim == 3
    set(gcf,'Position',[0 0 800 600], 'color', 'w');
    set(gca,'Position',[0.01 0.01 0.99,0.99]);
end
Str = {};
h = [];
if nDim == 3
    hold on;
    h = [h plot3(Ref0(1, :), Ref0(2, :), Ref0(3, :), 'g.')];
    Str{end+1} = 'target';
else
    plot(Ref0(1, :), Ref0(2, :), 'g.');
end
optR = TModel(end).Tf(1:nDim, 1:nDim);
optT = TModel(end).Tf(1:nDim, end);
% Aft = Loc2Glo(Mov0, optR', optT);
% [NNIdx, DD] = knnsearch(Ref0', Aft');
% tmp = Aft - Ref0(:, NNIdx);
% Err = sqrt(sum(tmp.^2));
Err = compute_residualFun(params); 
Z = expectationFun(Err, optModel);
if size(Z, 1) > 1
    [~, Label] = max(Z);
else
    Label = Z;
end
P = cat(2, optModel(:).prior)
Coeff = cat(2, optModel(:).coeff)
C = unique(Label);
ColorArray = {'r', 'b', 'k', 'm'};
for id = 1 : 1 : length(C)
    selId = C(id);
    idx = find(Label == selId);
    colorId = mod(C(id), length(ColorArray));
    if colorId == 0
        colorId = length(ColorArray);
    end
    data = Aft(:, idx);
    color = ColorArray{colorId};
    if nDim == 2
        ax = plot(data(1, :), data(2, :), 'color', color, 'marker', 'o', 'markersize', 3, 'linestyle', 'none' );
    else
        pcshow(data', color, 'markersize', 50);
    end
    str = sprintf('source with n_i = %.2f, coeff = %.2f', optModel(selId).prior, optModel(selId).coeff );
    Str{end+1} = str;
    hold on;
end
legend(Str, 'box', 'off', 'location', 'best', 'FontSize', 16);
xlabel('X/m', 'FontSize', 16);
ylabel('Y/m', 'FontSize', 16);
zlabel('Z/m', 'FontSize', 16);
title('Sparsity Analysis', 'FontSize', 16);

figure;
hold on;
grid on;
axis equal;
if nDim == 3
    set(gcf,'Position',[0 0 800 600], 'color', 'w');
    set(gca,'Position',[0.01 0.01 0.99,0.99]);
    view(3);
    showPointCloud(Ref0', 'g', 'markersize', 100 );
    showPointCloud(Aft', 'b',  'markersize', 100 );
    % showPointCloud(Mov0', 'r', 'markersize', 100 );
end
if nDim == 2
    plot(Ref0(1, :), Ref0(2, :), 'g.');
   % plot(Mov0(1, :), Mov0(2, :), 'r.');
    plot(Aft(1, :), Aft(2, :), 'bo', 'markersize', 3);
end
title('MiNoM Results');
% NoiseAnalysisFun(optModel, Err); 
end
