function varargout = MoEPFittingMultiFun(Err, P0, maxEMCounter, maxIter_em, bRand, IS_SHOW)
% [f_empirical, x_empirical] = ecdf(Err);
[N, edges] = histcounts(Err, 100);
gap = edges(2) - edges(1);
h_emprical = N' / sum(N);
cdf_empirical = cumsum(h_emprical);
% center = (edges(1:end-1) + edges(2:end)) / 2;
% figure;
% hold on;
% bar( center, h_emprical);
scale = 1000.0;
TWD = [];
TModel = [];
for id = 1 : 1 : maxEMCounter
    [Model, Z] = MoEPFittingFun(Err, P0, maxIter_em, bRand);
    % cat(2, Model(:).para)
    aa = cdfMoEPFun(Model, edges');
    h_theory = aa(2:end) - aa(1:end-1);
    cdf_theory    = cumsum(h_theory);
    h_diff = abs(cdf_theory - cdf_empirical);
    WD = scale * sum(h_diff) * gap;
    TWD(end+1) = WD;
    tmp = [];
    tmp.model = Model;
    tmp.z     = Z;
    TModel = [TModel tmp];
end
[~, id] = min(TWD);
varargout{1} = TModel(id).model;
if nargout == 2
    varargout{2} = TModel(id).z;
end
if IS_SHOW
    figure;
    hold on;
    grid on;
    plot(TWD, 'bo');
    [~, id] = min(TWD);
    NoiseAnalysisFun(TModel(id).model, Err);
    title('Optimal');
    [~, id] = max(TWD);
    NoiseAnalysisFun(TModel(id).model, Err);
    title('Worst');
end
end
