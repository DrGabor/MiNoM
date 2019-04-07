function [varargout] = compute_residual(params)
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
if strcmpi(params.mode, 'plane2plane')
    Sigma_Mov = params.Sigma_Mov; 
    Sigma_Ref = params.Sigma_Ref; 
    Tf = params.Tf; 
    dR = Tf(1:3, 1:3); 
    % dT = Tf(1:3, end); 
    MArray = [];
    parfor id = 1 : 1 : length(Aft)
        cov_mov = dR * Sigma_Mov(id).cov * dR';
        cov_ref = Sigma_Ref(id).cov;
        M = inv(cov_mov + cov_ref);
        m = [M(1, 1) M(1, 2) M(1, 3) M(2, 2) M(2, 3) M(3, 3)]';
        MArray = [MArray m];
        Err = [Err sqrt(V(:, id)'*M*V(:, id))];
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