function doses = getStats(wellMask, doseMap, includePlots)
% Asumimos que tienen las mismas dimensiones !(!!!!)
doseMap.resize([wellMask.minX, wellMask.maxX, wellMask.minY, wellMask.maxY]);

doses = doseMap.data(wellMask.data==1);
doseMapCopy = doseMap.copy;
doseMapCopy.data(wellMask.data==0) = nan;
meanDose = mean(doses);
stdDose = std(doses);

if exist('includePlots') && includePlots==true
    % Make plots
    figure(1);
    subplot(2,2,1);
    doseMapCopy.plotLinear('X',0);
    axis([-0.5 0.5 0 max(doseMap.getLinear('X',0))]);
    title('X');
    subplot(2,2,2);
    doseMapCopy.plotLinear('Y',0);
    axis([-0.5 0.5 0 max(doseMap.getLinear('Y',0))]);
    title('Y');
    subplot(2,1,2);
    doseMap.plotSlice;
    colorbar
    hold on
    contour(wellMask.getAxisValues('X'),wellMask.getAxisValues('Y'),wellMask.data);
    axis([-1 1 -1 1]);
    
    % Make DVH
    figure(2)
    dvh = nan(100, 1);
    level100 = meanDose;
    doseValues = 0:150;
    for i=1:numel(doseValues)
        dvh(i) = 100 * sum(doses > level100*(i/100)) / numel(doses);
    end
    plot(doseValues, dvh);
    xlabel('Relative dose');
    ylabel('Relative volume');
    axis([0 150 0 100]);
    
    % Make dose histogram
    figure(3);
    hist(doses,meanDose*(0.6:0.1:1.5));
    title('Relative dose in well');
    
    relDoses = 100 * doses / meanDose;
    r98_102 = sum(relDoses>=98 & relDoses<=102) / numel(relDoses);
    r95_105 = sum(relDoses>=95 & relDoses<=105) / numel(relDoses);
    r90_110 = sum(relDoses>=90 & relDoses<=110) / numel(relDoses);
    r80_120 = sum(relDoses>=80 & relDoses<=120) / numel(relDoses);
end
end

