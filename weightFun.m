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