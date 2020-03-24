clear all
load('allSigmas.mat');

%% Fit and plot

polyFlu = nan(maxE, 2);
polyEdep = nan(maxE, 2);

for i=minE:maxE
    
    allSteps = 1:NZ;
    modelFun = @(b,x) b(1)*x.^2 + b(2)*x;

    figure(i)
    subplot(2,1,1);
    
    minVal = find(~isnan(sigFlu{i}),1);
    maxVal = find(~isnan(sigFlu{i}),1,'last');
    
    % Probar el peso como (1 + x - minVal).^(-2);
    w = (1 + minVal + allSteps).^(-2);
    w(isnan(sigFlu{i})) = nan;
    % plot(w)
    
    % Fit
    X = (ZValues(minVal:maxVal))';
    Y = (sigFlu{i}(minVal:maxVal))';
    F = fit(X, Y, 'poly2', 'Lower', [0 0 0], 'Upper', [Inf Inf 0]);
    polyF = [F.p1 F.p2];
    nlm = fitnlm(X,Y,modelFun,polyF,'Weight',w(minVal:maxVal));
    polyW = nlm.Coefficients.Estimate';
    polyFlu(i, :) = polyW;
    
    % Plot
    plot(X,Y,'b.');
    hold on
    plot(ZValues(1:maxVal), modelFun(polyF,ZValues(1:maxVal)), 'r-');
    plot(ZValues(1:maxVal), modelFun(polyW,ZValues(1:maxVal)), 'g-');
    
    subplot(2,1,2);
    
    minVal = find(~isnan(sigEdep{i}),1);
    maxVal = find(~isnan(sigEdep{i}),1,'last');
    
    % Probar el peso como (1 + x - minVal).^(-2);
    w = (1 + minVal + allSteps).^(-2);
    w(isnan(sigEdep{i})) = nan;
    % plot(w)
    
    % Fit
    X = (ZValues(minVal:maxVal))';
    Y = (sigEdep{i}(minVal:maxVal))';
    F = fit(X, Y, 'poly2', 'Lower', [0 0 0], 'Upper', [Inf Inf 0]);
    polyF = [F.p1 F.p2];
    nlm = fitnlm(X,Y,modelFun,polyF,'Weight',w(minVal:maxVal));
    polyW = nlm.Coefficients.Estimate';
    polyEdep(i, :) = polyW;
    
    % Plot
    plot(X,Y,'b.');
    hold on
    plot(ZValues(1:maxVal), modelFun(polyF,ZValues(1:maxVal)), 'r-');
    plot(ZValues(1:maxVal), modelFun(polyW,ZValues(1:maxVal)), 'g-');    
end

save('polyFlu.mat','polyFlu');
save('polyEdep.mat','polyEdep');