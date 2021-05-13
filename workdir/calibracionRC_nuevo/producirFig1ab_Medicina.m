%% Abrir todos los archivos
clear all; close all

% Choose dataset
%filmIndex = 3; % EBT2,2=EBT3,3=EBT3unl
dataIndex = 3; % 1=HF fotones, 2=FQS protones, 3=Medicina
figure(1)
set(gcf, 'Position',  [100, 100, 900, 600])
ha = tight_subplot(2,2,0.02,[0.07 0.05],[0.06 0.01])

for filmIndex=1:2
    axes(ha(filmIndex))
    %% Load data
    if dataIndex==1
        load('basicData_HF.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
    elseif dataIndex==2
        load('basicData_FQS.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
    elseif dataIndex==3
        load('basicData_Medicina_SCANNED.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
    else
        error('No data available.');
    end
    
    %% Create single value table
    meanValues = nan(filmsPerSample, Nsamples, 3);
    stdValues = nan(filmsPerSample, Nsamples, 3);
    for i=1:filmsPerSample
        for j=1:Nsamples
            for l=1:3
                meanValues(i,j,l) = mean2(imageSubsets{i,j}(:,:,l));
                stdValues(i,j,l) = std2(imageSubsets{i,j}(:,:,l));
            end
        end
    end
    
    finalPixelValues = nan(filmsPerSample, Nsamples, 3);
    finalPixelErrValues = nan(filmsPerSample, Nsamples, 3);
    for i=1:filmsPerSample
        for j=1:Nsamples
            for l=1:3
                finalPixelValues(i,j,l) = mean(meanValues(i,j,l)) / maxBits;
                finalPixelErrValues(i,j,l) = rssq(stdValues(i,j,l)) / scansPerSample / maxBits;
            end
        end
    end
    
    %% Estudiar SOLAMENTE uno de los tipos
    Rmask = 1:10;
    Bmask = 1:10;
    Gmask = 1:10;
    
    R = finalPixelValues(filmIndex,:,1);
    G = finalPixelValues(filmIndex,:,2);
    B = finalPixelValues(filmIndex,:,3);
    dR = finalPixelErrValues(filmIndex,:,1);
    dG = finalPixelErrValues(filmIndex,:,2);
    dB = finalPixelErrValues(filmIndex,:,3);
    dosePoints = 0:0.1:(round(max(dosesGy)));
    
    %% Fit type I
    figure(1)
    
    [fR1, gofR1] = fitTypeI(dosesGy(Rmask), R(Rmask), dR(Rmask));
    plot(dosePoints, fR1(dosePoints), 'r-'); hold on
    [fG1, gofG1] = fitTypeI(dosesGy, G, dG);
    plot(dosePoints, fG1(dosePoints), 'g-');
    [fB1, gofB1] = fitTypeI(dosesGy(Bmask), B(Bmask), dB(Bmask));
    plot(dosePoints, fB1(dosePoints), 'b-');
    title(sprintf('%s', filmOrder{filmIndex}));
    errorbar(dosesGy, R, dR, 'r.', 'MarkerSize', 10)
    hold on
    errorbar(dosesGy, G, dG, 'g.', 'MarkerSize', 10)
    errorbar(dosesGy, B, dB, 'b.', 'MarkerSize', 10)
    grid on
    %xlabel('Doses (Gy)')
    if (filmIndex==1)
        ylabel('Relative PV');
    end
    xlim([0 round(max(dosesGy))]);
    
    CoefR1 = [fR1.alpha fR1.beta fR1.gamma];
    CoefG1 = [fG1.alpha fG1.beta fG1.gamma];
    CoefB1 = [fB1.alpha fB1.beta fB1.gamma];
    
    ypos = max([max(R), max(G), max(B)]) * 0.8;
    %text(8,ypos,sprintf('R^2 = %3.4f',gofR1.rsquare),'Color','red')
    %text(8,ypos*0.9,sprintf('R^2 = %3.4f',gofG1.rsquare),'Color','green')
    %text(8,ypos*0.8,sprintf('R^2 = %3.4f',gofB1.rsquare),'Color','blue')
    
    axis([0 15 0 1]);
    
    %% Fit Type 3
    NODR = log10(R(1) ./ R);
    dNODR = (1./R./log(10)).*dR;
    NODG = log10(G(1) ./ G);
    dNODG = (1./G./log(10)).*dG;
    NODB = log10(B(1) ./ B);
    dNODB = (1./B./log(10)).*dB;
    
    
    if (filmIndex==1)
        ylabel('netOD');
    end
    grid on
    
     axes(ha(filmIndex+2))
    
    [fR3, gofR3] = fitTypeIII_old(dosesGy(Rmask), NODR(Rmask), dNODR(Rmask));
    plot(dosePoints, fR3(dosePoints), 'r-'); hold on
    [fG3, gofG3] = fitTypeIII_old(dosesGy, NODG, dNODG);
    plot(dosePoints, fG3(dosePoints), 'g-');
    [fB3, gofB3] = fitTypeIII_old(dosesGy(Bmask), NODB(Bmask), dNODB(Bmask));
    plot(dosePoints, fB3(dosePoints), 'b-');
    
    errorbar(dosesGy, NODR, dNODR, 'r.', 'MarkerSize', 10)
    hold on
    errorbar(dosesGy, NODG, dNODG, 'g.', 'MarkerSize', 10)
    errorbar(dosesGy, NODB, dNODB, 'b.', 'MarkerSize', 10)
    
    %title(sprintf('%s', filmOrder{filmIndex}));
    %text(1,0.5,sprintf('R^2 = %3.4f',gofR3.rsquare),'Color','red')
    %text(1,0.45,sprintf('R^2 = %3.4f',gofG3.rsquare),'Color','green')
    %text(1,0.4,sprintf('R^2 = %3.4f',gofB3.rsquare),'Color','blue')
    ylim([0 1]);
    
    CoefR3 = [fR3.alpha fR3.beta fR3.gamma];
    CoefG3 = [fG3.alpha fG3.beta fG3.gamma];
    CoefB3 = [fB3.alpha fB3.beta fB3.gamma];
end
