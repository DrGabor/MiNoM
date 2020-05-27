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
