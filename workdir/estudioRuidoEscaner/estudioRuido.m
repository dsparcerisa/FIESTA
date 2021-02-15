%% Cargar imágenes
clear all; close all
basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_08_12 - Medida campo escaner nuevo/Modo VIDRIERA con cristales';
scansPerSample = 5;
%filmCode = {'Empty', 'EBT2', 'EBT3', 'EBT3unl'};
filmCode = {'Empty', 'EBT3b', 'EBT3unl'};
cropRect = [400 650 1000 1440];
for j=1:numel(filmCode)
    filmTypeIndex = j;
    
    figure(j)
    set(gcf, 'Position',[10 10 1100 400]);
    %% Read all images
    I = {};
    for i=1:scansPerSample
        fileName = sprintf('%s_%i.tif',filmCode{filmTypeIndex},i);
        filePath = [basePath filesep fileName];
        I{i} = imread(filePath);
        %I{i} = imcrop(I{i}, cropRect);
    end
    [pixelsXcm, maxBits] = getImgMetaInfo(filePath);
    dxy = 1/pixelsXcm;
    
    %% Merge all images
    Imerged = zeros(size(I{1}));
    for i=1:scansPerSample
        Imerged = Imerged + double(I{i});
    end
    Imerged = Imerged / scansPerSample;
    
    %% Analyze and show
    [NY, NX, colors] = size(Imerged);
    Xvalues = dxy*(-((NX-1)/2):((NX-1)/2));
    Yvalues = dxy*(-((NY-1)/2):((NY-1)/2));
    subplot(1,4,1);
    imshow(uint16(Imerged));
    title(filmCode{filmTypeIndex});
    
    subplot(2,4,2);
    meanX = squeeze(mean(Imerged,1));
    meanXR = meanX(:,1); 
    meanXG = meanX(:,2);
    meanXB = meanX(:,3);
 
    % Quick fix
    meanXR((end-2):end) = nan;
    meanXG((end-2):end) = nan;
    meanXB((end-2):end) = nan;   
    
    plot(Xvalues,meanXR,'r'); hold on;
    plot(Xvalues,meanXG,'g')
    plot(Xvalues,meanXB,'b')
    title('Transversal direction');
    xlabel('Distance from center (cm)');    
    ylabel('PV');
    xlim([Xvalues(1) Xvalues(end)]);
    
    subplot(2,4,6);
    meanY = squeeze(mean(Imerged,2));
    meanYR = meanY(:,1);
    meanYG = meanY(:,2);
    meanYB = meanY(:,3);
    plot(Yvalues,meanYR,'r'); hold on;
    plot(Yvalues,meanYG,'g')
    plot(Yvalues,meanYB,'b')
    ylabel('PV');
    title('Scanning direction');
    xlabel('Distance from center (cm)');    
    xlim([Yvalues(1) Yvalues(end)]);
    
    subplot(2,4,3);
    plot(Xvalues,meanXR ./ max(meanXR),'r'); hold on;
    plot(Xvalues,meanXG ./ max(meanXG),'g')
    plot(Xvalues,meanXB ./ max(meanXB),'b')
    title('Transversal direction');
    xlabel('Distance from center (cm)');    
    ylabel('Normalized PV');
    xlim([Xvalues(1) Xvalues(end)]);
    ylim([0.9 1]);
    
    subplot(2,4,7);
    plot(Yvalues,meanYR ./ max(meanYR),'r'); hold on;
    plot(Yvalues,meanYG ./ max(meanYG),'g')
    plot(Yvalues,meanYB ./ max(meanYB),'b')
    ylabel('Normalized PV');
    title('Scanning direction');
    xlabel('Distance from center (cm)');    
    xlim([Yvalues(1) Yvalues(end)]);
    ylim([0.9 1]);    

    subplot(2,4,4);
    plot(Xvalues,meanXR - max(meanXR),'r'); hold on;
    plot(Xvalues,meanXG - max(meanXG),'g')
    plot(Xvalues,meanXB - max(meanXB),'b')
    title('Transversal direction');
    ylabel('PV - PV_m_a_x');
    xlabel('Distance from center (cm)');
    xlim([Xvalues(1) Xvalues(end)]);
    
    subplot(2,4,8);
    plot(Yvalues,meanYR - max(meanYR),'r'); hold on;
    plot(Yvalues,meanYG - max(meanYG),'g')
    plot(Yvalues,meanYB - max(meanYB),'b')
    title('Scanning direction');
    ylabel('PV - PV_m_a_x');
    xlabel('Distance from center (cm)');    
    xlim([Yvalues(1) Yvalues(end)]);
    
end