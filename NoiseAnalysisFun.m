function [] = NoiseAnalysisFun(Model, Err, Edges )
if ~isrow(Err)
    Err = Err';
end
Tau = cat(2, Model(:).para);
Coeff = cat(2, Model.coeff);
P = cat(2, Model.prior);
hhh = histcounts(Err, Edges);
hhh = hhh / sum(hhh); 
%%
gtH = zeros(1, length(Edges));
for k = 1 : 1 : length(Coeff)
    tau = Tau(k);
    p   = P(k);
    Q   = gammainc(Edges.^p*tau, ones(1, length(Edges))/p);
    gtH = gtH + Coeff(k) * Q;
end
markerIndex = 1 : 10 : length(Edges); 
gtH = gtH(2:end) - gtH(1:end-1);
figure; 
hold on; 
grid on; 
box on; 
Center = ( Edges(1:1:end-1) + Edges(2:1:end))/2;
bar(Center,  hhh, 'r');
plot(Center, gtH, 'gp-', 'linewidth', 3, 'markersize', 5, 'markerindices', markerIndex);
set(gca, 'XLim', [min(Edges) max(Edges)]); 

legend({'Empirical Histogram', 'Theoretical Histogram'}, 'FontSize', 22, ...
    'FontWeight', 'bold', 'box', 'off', 'location', 'best' );
xlabel('Residual Error/m', 'FontSize', 22, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 22, 'FontWeight', 'bold'); 
title('RED Analysis', 'FontSize', 22);

end

