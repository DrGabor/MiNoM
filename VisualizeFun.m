function VisualizeFun(TModel, params)
Mov0 = params.Mov0; 
Ref0 = params.Ref0; 
Aft = params.Aft;
optModel = TModel(end).model;
tmp = struct('data', []);
P = cat(1, optModel(:).prior);
TSita = repmat(tmp, 1, length(P));
for id = 1 : 1 : length(TModel)
    optModel = TModel(id).model;
    for k = 1 : 1 : length(optModel)
        if strcmp(optModel(k).type, 'ep')
            p = optModel(k).prior;
            [~, ind] = ismember(p, P);
            TSita(ind).data(end+1) = optModel(k).para;
        end
    end
end
h = figure;
subplot(121);
hold on;
grid on;
llh = cat(1, TModel(:).llh);
plot(llh, 'b.-');
title('Likelihood Function');
xlabel('Iteration');
ylabel('llh');
%%%%%%%%%%
figure(h);
subplot(122);
hold on;
grid on;
Str = {};
h = [];
ColorArray = {'r', 'b', 'm', 'k', 'g', 'c'};
for id = 1 : 1 : length(TSita)
    y = TSita(id).data;
    colorId = mod(id, length(ColorArray));
    if colorId == 0
        colorId = length(ColorArray);
    end
    color = ColorArray{colorId};
    h = [h plot(y, 'marker', '.', 'linestyle', '-', 'color', color)];
    Str{end+1} = sprintf('p = %.1f', P(id));
end
xlabel('Iteration');
ylabel('$\bf{\theta}$ value', 'interpreter', 'latex');
legend(h, Str);
legend('boxoff');
legend('location', 'best');
str = 'covergence curve for $\bf{\theta}$';
title(str, 'interpreter', 'latex');
%%%%%%%%%%
%% plot registration results.
nDim = size(Ref0, 1);
figure;
hold on;
grid on;
% axis equal;
if nDim == 3
    set(gcf,'Position',[0 0 800 600], 'color', 'w');
    set(gca,'Position',[0.01 0.01 0.99,0.99]);
end
Str = {};
h = [];
if nDim == 3
    hold on;
    h = [h plot3(Ref0(1, :), Ref0(2, :), Ref0(3, :), 'g.')];
    Str{end+1} = 'target';
else
    plot(Ref0(1, :), Ref0(2, :), 'g.');
end
optR = TModel(end).Tf(1:nDim, 1:nDim);
optT = TModel(end).Tf(1:nDim, end);
% Aft = Loc2Glo(Mov0, optR', optT);
% [NNIdx, DD] = knnsearch(Ref0', Aft');
% tmp = Aft - Ref0(:, NNIdx);
% Err = sqrt(sum(tmp.^2));
Err = compute_residualFun(params); 
Z = expectationFun(Err, optModel);
if size(Z, 1) > 1
    [~, Label] = max(Z);
else
    Label = Z;
end
P = cat(2, optModel(:).prior)
Coeff = cat(2, optModel(:).coeff)
C = unique(Label);
ColorArray = {'r', 'b', 'k', 'm'};
for id = 1 : 1 : length(C)
    selId = C(id);
    idx = find(Label == selId);
    colorId = mod(C(id), length(ColorArray));
    if colorId == 0
        colorId = length(ColorArray);
    end
    data = Aft(:, idx);
    color = ColorArray{colorId};
    if nDim == 2
        ax = plot(data(1, :), data(2, :), 'color', color, 'marker', 'o', 'markersize', 3, 'linestyle', 'none' );
    else
        pcshow(data', color, 'markersize', 50);
    end
    str = sprintf('source with n_i = %.2f, coeff = %.2f', optModel(selId).prior, optModel(selId).coeff );
    Str{end+1} = str;
    hold on;
end
legend(Str, 'box', 'off', 'location', 'best', 'FontSize', 16);
xlabel('X/m', 'FontSize', 16);
ylabel('Y/m', 'FontSize', 16);
zlabel('Z/m', 'FontSize', 16);
title('Sparsity Analysis', 'FontSize', 16);

figure;
hold on;
grid on;
axis equal;
if nDim == 3
    set(gcf,'Position',[0 0 800 600], 'color', 'w');
    set(gca,'Position',[0.01 0.01 0.99,0.99]);
    view(3);
    showPointCloud(Ref0', 'g', 'markersize', 100 );
    showPointCloud(Aft', 'b',  'markersize', 100 );
    % showPointCloud(Mov0', 'r', 'markersize', 100 );
end
if nDim == 2
    plot(Ref0(1, :), Ref0(2, :), 'g.');
   % plot(Mov0(1, :), Mov0(2, :), 'r.');
    plot(Aft(1, :), Aft(2, :), 'bo', 'markersize', 3);
end
title('MiNoM Results');
% NoiseAnalysisFun(optModel, Err); 
end