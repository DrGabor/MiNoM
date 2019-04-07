function [f1] = cdfMoEPFun(Model,xx)
f1 = zeros(length(xx), 1);
for id = 1 : 1 : length(Model)
    p = Model(id).prior;
    sita = Model(id).para;
    w = Model(id).coeff;
    f1 = f1 + w * gammainc(xx.^p*sita, ones(length(xx), 1)/p);
end
end

