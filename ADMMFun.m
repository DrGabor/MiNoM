function [ varargout ] = ADMM( params, Z0, Model )
Z = [];
Sita = [];
P = [];
for k = 1 : 1 : length(Model)
    if strcmp(Model(k).type, 'ep')
        Z = [Z; Z0(k, :)];
        Sita = [Sita; Model(k).para];
        P = [P; Model(k).prior];
    end
end
par = [];
par.mu = 10.0;
par.maxMu = 1e5;
par.maxInner = 5;
par.maxOuter = 20;
%%%%%%%%%% Using Lagrange Multipliers to solve the Lp Norm optimization.
miu = par.mu;
W_P = Z .* repmat(Sita, 1, length(Z));
Aft = params.Aft;
Mov = params.Mov;
Ref = params.Ref;
Ref_Normal = params.Ref_Normal;
if strcmpi(params.mode, 'point2plane')
    Lamda = zeros(1, size(Ref, 2) );
else
    Lamda = zeros(size(Ref));
end
R = eye(3);
T = zeros(3, 1);
if strcmpi(params.mode, 'plane2plane')
    MArray = params.MArray;
    ss = struct('data', []);
    LArray = repmat(ss, 1, length(MArray));
    parfor id = 1 : 1 : length(MArray)
        m = MArray(:, id);
        %m = [M(1, 1) M(1, 2) M(1, 3) M(2, 2) M(2, 3) M(3, 3)]';
        M = [m(1) m(2) m(3)
            m(2) m(4) m(5)
            m(3) m(5) m(6)];
        L = chol(M);
        invL = inv(L);
        LArray(id).L = L;
        LArray(id).invL = invL;
    end
end
for outer = 1 : 1 : par.maxOuter
    dual = 0.0;
    %%%%%%%%%%%%% this loop instaneously recovery the rotation and
    %%%%%%%%%%%%% translate.
    for inner = 1 : 1 : par.maxInner
        if strcmpi(params.mode, 'point2point')
            U = Aft - Ref;
        end
        if strcmpi(params.mode, 'point2plane')
            U =  sum(Ref_Normal .* (Aft - Ref)); 
        end
        if strcmpi(params.mode, 'plane2plane')
            Diff = Aft - Ref;
            for id = 1 : 1 : length(LArray)
                U(:, id) = LArray(id).L * Diff(:, id);
            end
        end
        % step2.1 update residual error vector.
        H = U + Lamda / miu;
        selId = find(P == 2);
        W2 = W_P(selId, :);
        Factor = miu ./ (miu + 2*W2);
        selId = find(P == 1);
        W1 = W_P(selId, :);
        V = Factor .* H;
        VNorm = VectorNorm(V);
        Factor = W1 ./ (miu + 2*W2);
        UnitV = V ./ VNorm;
        Factor = max([VNorm - Factor; zeros(1, length(VNorm))]);
        E = Factor .* UnitV; % zeros(size(U));
        % step2.2 update transformation.
        params.Aft = Aft;
        params.W = ones(1, length(Aft));
        if strcmpi(params.mode, 'point2point')
            Ref_Modify = Ref + E - Lamda / miu;
            params.Ref = Ref_Modify;
            [dR dT] = reg_point2point(params);
        end
        if strcmpi(params.mode, 'point2plane')
            Ref_Modify = Ref + (E - Lamda / miu) .* Ref_Normal;
            params.Ref = Ref_Modify;
            [dR dT] = reg_point2plane(params);
        end
        
        if strcmpi(params.mode, 'plane2plane')
            Ref_Modify = zeros(size(Ref));
            for id = 1 : 1 : length(Ref_Modify)
                Ref_Modify(:, id) = Ref(:, id) + LArray(id).invL * (E(:, id) - Lamda(:, id) / miu);
            end
            params.Ref = Ref_Modify;
            [dR dT] = reg_plane2plane(params);
        end
        Aft_New = Loc2Glo(Aft, dR', dT);
        tmpNorm = VectorNorm(Aft - Aft_New);
        dual = max(tmpNorm);
        Aft = Aft_New;
        %%%%%%%%%%%%% move dr and dt to global coordinate.
        R = dR * R;
        T = dR * T + dT;
        if dual < 1e-4
            break;
        end
    end
    % step2.3 update Lagrange multipliers.
    Lamda = Lamda + miu * (U - E);
    if miu < par.maxMu
        miu = 1.2 * miu;
    end
    primal = max(VectorNorm(U - E));
    if primal < 1e-4 & dual < 1e-4
        break;
    end
end
if nargout == 1
    varargout{1} =  XArray;
end
if nargout == 2
    varargout{1} = R;
    varargout{2} = T;
end
bTest = 1;
end