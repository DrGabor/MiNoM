function Model = maximation( Err, Z, Model0 )
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