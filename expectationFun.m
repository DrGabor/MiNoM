function Z = expectation( Err, Model )
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