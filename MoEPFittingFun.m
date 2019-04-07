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
