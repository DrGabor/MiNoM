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

