function [ varargout ] = IRLS( params, Z0, Model, maxIter )
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
bTest = 1;
end