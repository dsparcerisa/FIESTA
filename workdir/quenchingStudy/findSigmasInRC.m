function [finalSigmas,dfinalSigmas] = findSigmasInRC(allI,NValidPoints,CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits, directTo)
finalSigmas = nan(NValidPoints,1);
dfinalSigmas = nan(NValidPoints,1);

if nargin<9
    directTo = 1:NValidPoints;
end

for i=directTo
    d1 = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    d2 = getDoseMicke(double(allI{i+11}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    
    % Probar con la suma
    sumX = mean(d1.data,2);
    baseLineX = min(sumX);
    sumX = sumX - baseLineX;
    X = d1.getAxisValues('X');
    subplot(4,2,1); hold off
    plot(X, sumX, 'b.')
    hold on
    title('X low (1)');
    [FX, gofx] = fit(X', sumX, 'gauss1');
    SX = 10*FX.c1/sqrt(2);
    FXerrors = confint(FX);
    dSX = 10*max(abs(FXerrors(:,3)-FX.c1)) / sqrt(2);
    plot(FX)
    
    sumY = mean(d1.data);
    baseLineY = min(sumY);
    sumY = sumY - baseLineY;
    Y = d1.getAxisValues('Y');
    subplot(4,2,2); hold off
    plot(Y, sumY, 'b.')
    hold on
    title('Y low (2)');
    [FY, gofy] = fit(Y', sumY', 'gauss1');
    SY = 10*FY.c1/sqrt(2);
    FYerrors = confint(FY);
    dSY = 10*max(abs(FYerrors(:,3)-FY.c1)) / sqrt(2);
    plot(FY)
    
    sumX2 = mean(d2.data,2);
    baseLineX2 = min(sumX2);
    sumX2 = sumX2 - baseLineX2;
    X2 = d2.getAxisValues('X');
    subplot(4,2,3); hold off
    plot(X2, sumX2, 'b.')
    hold on
    title('X high (4)');
    [FX2, gofx2] = fit(X2', sumX2, 'gauss1');
    SX2 = 10*FX2.c1/sqrt(2);
    FX2errors = confint(FX2);
    dSX2 = 10*max(abs(FX2errors(:,3)-FX2.c1)) / sqrt(2);
    plot(FX2)
    
    sumY2 = mean(d2.data);
    baseLineY2 = min(sumY2);
    sumY2 = sumY2 - baseLineY2;
    Y2 = d2.getAxisValues('Y');
    subplot(4,2,4); hold off
    plot(Y2, sumY2, 'b.')
    hold on
    title('Y high (5)');
    [FY2, gofy2] = fit(Y2', sumY2', 'gauss1');
    SY2 = 10*FY2.c1/sqrt(2);
    FY2errors = confint(FY2);
    dSY2 = 10*max(abs(FY2errors(:,3)-FY2.c1)) / sqrt(2);
    plot(FY2)
    
    SXY = (SX / dSX / dSX + SY / dSY / dSY) / (1 / dSX / dSX + 1 / dSY / dSY);
    dSXY = (1 / dSX / dSX + 1 / dSY / dSY)^(-0.5);
    
    SXY2 = (SX2 / dSX2 / dSX2 + SY2 / dSY2 / dSY2) / (1 / dSX2 / dSX2 + 1 / dSY2 / dSY2);
    dSXY2 = (1 / dSX2 / dSX2 + 1 / dSY2 / dSY2)^(-0.5);
    
    % Doble gaussiana
    [FX_DG, gofx_DG] = fit(X', sumX, 'gauss2', 'Lower', [0 -Inf 0 0 -Inf 0])
    FXerrors_DG = confint(FX_DG);   
    if FX_DG.a1 > FX_DG.a2
        SX_DG = 10*FX_DG.c1/sqrt(2);
        dSX_DG = 10*max(abs(FXerrors_DG(:,3)-FX_DG.c1)) / sqrt(2);
    else
        SX_DG = 10*FX_DG.c2/sqrt(2);
        dSX_DG = 10*max(abs(FXerrors_DG(:,6)-FX_DG.c2)) / sqrt(2);
    end
    subplot(4,2,5); hold off
    plot(X, sumX, 'b.');
    hold on
    plot(FX_DG)
    title('X low DG (7)');
    
    [FY_DG, gofy_DG] = fit(Y', sumY', 'gauss2', 'Lower', [0 -Inf 0 0 -Inf 0])
    FYerrors_DG = confint(FY_DG);   
    if FY_DG.a1 > FY_DG.a2
        SY_DG = 10*FY_DG.c1/sqrt(2);
        dSY_DG = 10*max(abs(FYerrors_DG(:,3)-FY_DG.c1)) / sqrt(2);
    else
        SY_DG = 10*FY_DG.c2/sqrt(2);
        dSY_DG = 10*max(abs(FYerrors_DG(:,6)-FY_DG.c2)) / sqrt(2);
    end
    subplot(4,2,6); hold off
    plot(Y, sumY, 'b.');
    hold on
    plot(FY_DG)
    title('Y low DG (8)');    
    
    [FX_DG2, gofx_DG2] = fit(X2', sumX2, 'gauss2', 'Lower', [0 -Inf 0 0 -Inf 0])
    FXerrors_DG2 = confint(FX_DG2);   
    if FX_DG2.a1 > FX_DG2.a2
        SX_DG2 = 10*FX_DG2.c1/sqrt(2);
        dSX_DG2 = 10*max(abs(FXerrors_DG2(:,3)-FX_DG2.c1)) / sqrt(2);
    else
        SX_DG2 = 10*FX_DG2.c2/sqrt(2);
        dSX_DG2 = 10*max(abs(FXerrors_DG2(:,6)-FX_DG2.c2)) / sqrt(2);
    end
    subplot(4,2,7); hold off
    plot(X2, sumX2, 'b.');
    hold on
    plot(FX_DG2)
    title('X high DG (10)');
    
    [FY_DG2, gofy_DG2] = fit(Y2', sumY2', 'gauss2', 'Lower', [0 -Inf 0 0 -Inf 0])
    FYerrors_DG2 = confint(FY_DG2);   
    if FY_DG2.a1 > FY_DG2.a2
        SY_DG2 = 10*FY_DG2.c1/sqrt(2);
        dSY_DG2 = 10*max(abs(FYerrors_DG2(:,3)-FY_DG2.c1)) / sqrt(2);
    else
        SY_DG2 = 10*FY_DG2.c2/sqrt(2);
        dSY_DG2 = 10*max(abs(FYerrors_DG2(:,6)-FY_DG2.c2)) / sqrt(2);
    end
    subplot(4,2,8); hold off
    plot(Y2, sumY2, 'b.');
    hold on
    plot(FY_DG2)
    title('Y high DG (11)');      

    SXY_DG = (SX_DG / dSX_DG / dSX_DG + SY_DG / dSY_DG / dSY_DG) / (1 / dSX_DG / dSX_DG + 1 / dSY_DG / dSY_DG);
    dSXY_DG = (1 / dSX_DG / dSX_DG + 1 / dSY_DG / dSY_DG)^(-0.5);
    
    SXY_DG2 = (SX_DG2 / dSX_DG2 / dSX_DG2 + SY_DG2 / dSY_DG2 / dSY_DG2) / (1 / dSX_DG2 / dSX_DG2 + 1 / dSY_DG2 / dSY_DG2);
    dSXY_DG2 = (1 / dSX_DG2 / dSX_DG2 + 1 / dSY_DG2 / dSY_DG2)^(-0.5);
    
    fprintf('\tSigma\tdSigma\tR2\n');
    fprintf('1 Xlo\t%3.4f\t%3.4f\t%3.4f\n', SX, dSX, gofx.rsquare);
    fprintf('2 Ylo\t%3.4f\t%3.4f\t%3.4f\n', SY, dSY, gofy.rsquare);
    fprintf('3 XYlo\t%3.4f\t%3.4f\t%3.4f\n', SXY, dSXY, 0);
    fprintf('4 Xhi\t%3.4f\t%3.4f\t%3.4f\n', SX2, dSX2, gofx2.rsquare);
    fprintf('5 Yhi\t%3.4f\t%3.4f\t%3.4f\n', SY2, dSY2, gofy2.rsquare);
    fprintf('6 XYhi\t%3.4f\t%3.4f\t%3.4f\n', SXY2, dSXY2, 0);
    fprintf('7 Xlo\t%3.4f\t%3.4f\t%3.4f\n', SX_DG, dSX_DG, gofx_DG.rsquare);
    fprintf('8 Ylo\t%3.4f\t%3.4f\t%3.4f\n', SY_DG, dSY_DG, gofy_DG.rsquare);
    fprintf('9 XYlo\t%3.4f\t%3.4f\t%3.4f\n', SXY_DG, dSXY_DG, 0);
    fprintf('10 Xhi\t%3.4f\t%3.4f\t%3.4f\n', SX_DG2, dSX_DG2, gofx_DG2.rsquare);
    fprintf('11 Yhi\t%3.4f\t%3.4f\t%3.4f\n', SY_DG2, dSY_DG2, gofy_DG2.rsquare);
    fprintf('12 XYhi\t%3.4f\t%3.4f\t%3.4f\n', SXY_DG2, dSXY_DG2, 0);    
    
    input_hf = input('Which one do we use: ');
    switch(input_hf)
        case 1
            finalSigmas(i) = SX;
            dfinalSigmas(i) = dSX;
        case 2
            finalSigmas(i) = SY;
            dfinalSigmas(i) = dSY;
        case 3
            finalSigmas(i) = SXY;
            dfinalSigmas(i) = dSXY;
        case 4
            finalSigmas(i) = SX2;
            dfinalSigmas(i) = dSX2;
        case 5
            finalSigmas(i) = SY2;
            dfinalSigmas(i) = dSY2;
        case 6
            finalSigmas(i) = SXY2;
            dfinalSigmas(i) = dSXY2;
        case 7
            finalSigmas(i) = SX_DG;
            dfinalSigmas(i) = dSX_DG;
        case 8
            finalSigmas(i) = SY_DG;
            dfinalSigmas(i) = dSY_DG;
        case 9
            finalSigmas(i) = SXY_DG;
            dfinalSigmas(i) = dSXY_DG;
        case 10
            finalSigmas(i) = SX_DG2;
            dfinalSigmas(i) = dSX_DG2;
        case 11
            finalSigmas(i) = SY_DG2;
            dfinalSigmas(i) = dSY_DG2;
        case 12
            finalSigmas(i) = SXY_DG2;
            dfinalSigmas(i) = dSXY_DG2;    
        otherwise
            error('Re-run this block');
    end
end
end

