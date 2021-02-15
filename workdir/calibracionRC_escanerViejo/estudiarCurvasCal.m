%% Abrir todos los archivos
clear all; close all
load('basicData.mat');

%% Create single value table
meanValues = nan(Nsamples, filmsPerSample, scansPerSample, 3);
stdValues = nan(Nsamples, filmsPerSample, scansPerSample, 3);
for i=1:Nsamples
   for j=1:filmsPerSample
       for k=1:scansPerSample
            for l=1:3
                meanValues(i,j,k,l) = mean2(imageSubsets{i,j,k}(:,:,l));
                stdValues(i,j,k,l) = std2(imageSubsets{i,j,k}(:,:,l));
            end
       end
   end
end

finalPixelValues = nan(Nsamples, filmsPerSample, 3);
finalPixelErrValues = nan(Nsamples, filmsPerSample, 3);
for i=1:Nsamples
    for j=1:filmsPerSample
        for l=1:3
            finalPixelValues(i,j,l) = mean(meanValues(i,j,:,l)) / maxBits;
            finalPixelErrValues(i,j,l) = rssq(stdValues(i,j,:,l)) / scansPerSample / maxBits;
        end
    end
end

%% Estudiar SOLAMENTE los EBT3-unl
filmIndex = 3;

Rmask = 1:10;
Bmask = 6:10;

R = finalPixelValues(:,filmIndex,1);
G = finalPixelValues(:,filmIndex,2);
B = finalPixelValues(:,filmIndex,3);
dR = finalPixelErrValues(:,filmIndex,1);
dG = finalPixelErrValues(:,filmIndex,2);
dB = finalPixelErrValues(:,filmIndex,3);
dosePoints = 0:0.1:12;

%% Fit type I
figure(1)
errorbar(dosesGy, R, dR, 'r.')
hold on
errorbar(dosesGy, G, dG, 'g.')
errorbar(dosesGy, B, dB, 'b.')
[fR1, gofR1] = fitTypeI(dosesGy(Rmask), R(Rmask), dR(Rmask));
plot(dosePoints, fR1(dosePoints), 'r-');
[fG1, gofG1] = fitTypeI(dosesGy, G, dG);
plot(dosePoints, fG1(dosePoints), 'g-');
[fB1, gofB1] = fitTypeI(dosesGy(Bmask), B(Bmask), dB(Bmask));
plot(dosePoints, fB1(dosePoints), 'b-');
[fB1b, gofB1b] = fit(dosesGy(Bmask)', B(Bmask), 'rat12', 'Weights', dB(Bmask).^(-2));
plot(dosePoints, fB1b(dosePoints), 'b:');
title(sprintf('Type I - Rational (%s)', filmOrder{filmIndex}));
grid on
xlabel('Doses (Gy)')
ylabel('Relative pixel values');
xlim([0 12]);

CoefR = [fR1.alpha fR1.beta fR1.gamma];
CoefG = [fG1.alpha fG1.beta fG1.gamma];
CoefB = [fB1.alpha fB1.beta fB1.gamma];

ypos = max([max(R), max(G), max(B)]) * 0.8;
text(8,ypos,sprintf('R^2 = %3.4f',gofR1.rsquare),'Color','red')
text(8,ypos*0.9,sprintf('R^2 = %3.4f',gofG1.rsquare),'Color','green')
text(8,ypos*0.8,sprintf('R^2 = %3.4f',gofB1.rsquare),'Color','blue')
text(8,ypos*0.7,sprintf('R^2 = %3.4f',gofB1b.rsquare),'Color','blue')

%% Fit type II
figure(2)

ODR = log10(maxBits ./ R);
dODR = (1./R./log(10)).*dR;
ODG = log10(maxBits ./ G);
dODG = (1./G./log(10)).*dG;
ODB = log10(maxBits ./ B);
dODB = (1./B./log(10)).*dB;

subplot(1,2,1);
errorbar(dosesGy, ODR, dODR, 'r.')
hold on
errorbar(dosesGy, ODG, dODG, 'g.')
errorbar(dosesGy, ODB, dODB, 'b.')
xlabel('Doses (Gy)')
ylabel('OD');
grid on
[fR2, gofR2] = fitTypeII(dosesGy(Rmask), ODR(Rmask), dODR(Rmask));
plot(dosePoints, fR2(dosePoints), 'r-');
[fG2, gofG2] = fitTypeII(dosesGy, ODG, dODG);
plot(dosePoints, fG2(dosePoints), 'g-');
[fB2, gofB2] = fitTypeII(dosesGy(Bmask), ODB(Bmask), dODB(Bmask));
plot(dosePoints, fB2(dosePoints), 'b-');

text(1,3.7,sprintf('R^2 = %3.4f',gofR2.rsquare),'Color','red')
text(1,3.6,sprintf('R^2 = %3.4f',gofG2.rsquare),'Color','green')
text(1,3.5,sprintf('R^2 = %3.4f',gofB2.rsquare),'Color','blue')

xlim([0 12]);
title(sprintf('Type II - Poly (%s)', filmOrder{filmIndex}));

subplot(1,2,2);
errorbar(dosesGy, ODR, dODR, 'r.')
hold on
errorbar(dosesGy, ODG, dODG, 'g.')
errorbar(dosesGy, ODB, dODB, 'b.')
xlabel('Doses (Gy)')
ylabel('OD');
grid on
[fR2r, gofR2r] = fitTypeIIb(dosesGy(Rmask), ODR(Rmask), dODR(Rmask));
plot(dosePoints, fR2r(dosePoints), 'r-');
[fG2r, gofG2r] = fitTypeIIb(dosesGy, ODG, dODG);
plot(dosePoints, fG2r(dosePoints), 'g-');
[fB2r, gofB2r] = fitTypeIIb(dosesGy(Bmask), ODB(Bmask), dODB(Bmask));
plot(dosePoints, fB2r(dosePoints), 'b-');

text(1,3.7,sprintf('R^2 = %3.4f',gofR2r.rsquare),'Color','red')
text(1,3.6,sprintf('R^2 = %3.4f',gofG2r.rsquare),'Color','green')
text(1,3.5,sprintf('R^2 = %3.4f',gofB2r.rsquare),'Color','blue')

xlim([0 12]);
title(sprintf('Type II - Rational (%s)', filmOrder{filmIndex}));

%% Type 3
NODR = log10(R(1) ./ R);
dNODR = (1./R./log(10)).*dR;
NODG = log10(G(1) ./ G);
dNODG = (1./G./log(10)).*dG;
NODB = log10(B(1) ./ B);
dNODB = (1./B./log(10)).*dB;

figure(3)
errorbar(dosesGy, NODR, dNODR, 'r.')
hold on
errorbar(dosesGy, NODG, dNODG, 'g.')
errorbar(dosesGy, NODB, dNODB, 'b.')
xlabel('Doses (Gy)')
ylabel('NOD');
grid on

[fR3, gofR3] = fitTypeIII(dosesGy(Rmask), NODR(Rmask), dNODR(Rmask));
plot(dosePoints, fR3(dosePoints), 'r-');
[fG3, gofG3] = fitTypeIII(dosesGy, NODG, dNODG);
plot(dosePoints, fG3(dosePoints), 'g-');
[fB3, gofB3] = fitTypeIII(dosesGy(Bmask), NODB(Bmask), dNODB(Bmask));
plot(dosePoints, fB3(dosePoints), 'b-');
title(sprintf('Type III - Potential (%s)', filmOrder{filmIndex}));text(1,1.1,sprintf('R^2 = %3.4f',gofR3.rsquare),'Color','red')
text(1,1.0,sprintf('R^2 = %3.4f',gofG3.rsquare),'Color','green')
text(1,0.9,sprintf('R^2 = %3.4f',gofB3.rsquare),'Color','blue')

%% Estudio de la solución de rat12
%syms pv p1 D p2 q1 q2
%rat12 = pv == (p1*D + p2)/(D^2 + q1*D + q2)
%S = solve(rat12, D)

fR1inv = @(pv) fR1.gamma + fR1.beta./(pv - fR1.alpha);
fG1inv = @(pv) fG1.gamma + fG1.beta./(pv - fG1.alpha);
fB1inv = @(pv) fB1.gamma + fB1.beta./(pv - fB1.alpha);

fB1binv_aux = @(pv,p1,p2,q1,q2) (p1 + (p1.^2 - 2.*p1.*pv.*q1 + pv.^2.*q1.^2 - 4.*q2.*pv.^2 + 4.*p2.*pv).^(1/2) - pv.*q1)./(2.*pv)
fB1binv = @(pv) real(fB1binv_aux(pv, fB1b.p1, fB1b.p2, fB1b.q1, fB1b.q2))
real([dosesGy' fR1inv(R) fG1inv(G) fB1inv(B) fB1binv(B)])

%% Estudio de autoconsistencia para todos:
figure(4)
dosesR = fR1inv(R);
ddosesR = max(abs(fR1inv(R+dR)-dosesR), abs(fR1inv(R-dR)-dosesR));
dosesG = fG1inv(G);
ddosesG = max(abs(fG1inv(G+dG)-dosesG), abs(fG1inv(G-dG)-dosesG));
dosesB = fB1inv(B);
ddosesB = max(abs(fB1inv(B+dB)-dosesB), abs(fB1inv(B-dB)-dosesB));
dosesB2 = fB1binv(B);
ddosesB2 = max(abs(fB1binv(B+dB)-dosesB2), abs(fB1binv(B-dB)-dosesB2));

plot(dosesGy, dosesGy, 'k-'); hold on

errorbar(dosesGy, dosesR, ddosesR, 'ro'); hold on
errorbar(dosesGy, dosesG, ddosesG, 'gx');
errorbar(dosesGy, dosesB, ddosesB, 'b+');
%errorbar(dosesGy, dosesB2, ddosesB2, 'bs');
%legend({'Red channel', 'Green channel', 'Blue channel-rat11', 'Blue channel-rat12'}, 'Location', 'SouthEast');
legend({'Reference', 'Red channel', 'Green channel', 'Blue channel'}, 'Location', 'SouthEast');

grid on
xlabel('Real dose (Gy)');
ylabel('Predicted dose (Gy)');
title(sprintf('Type I - Rational (%s)', filmOrder{filmIndex}));
resR = rssq((dosesGy(1:9)'-fR1inv(R(1:9))), 1) ./ 9;
resG = rssq((dosesGy'-fG1inv(G)), 1) / 10;
resB = rssq((dosesGy(Bmask)'-fB1inv(B(Bmask))), 1) / 6;
resB2 = rssq((dosesGy(Bmask)'-fB1binv(B(Bmask))), 1) / 6;
axis([0 12 0 12]);

text(2,10,sprintf('res = %3.2f Gy',resR),'Color','red', 'FontSize', 14)
text(2,9,sprintf('res = %3.2f Gy',resG),'Color','green', 'FontSize', 14)
text(2,8,sprintf('res = %3.2f Gy',resB),'Color','blue', 'FontSize', 14)
%text(2,7,sprintf('res = %3.2f Gy',resB2),'Color','blue', 'FontSize', 14)
set(gca, 'FontSize', 14)

%% Estudio con el Mayer para Tipo I
aR1 = @(pv) -fR1.beta ./ (pv - fR1.alpha).^2;
aG1 = @(pv) -fG1.beta ./ (pv - fG1.alpha).^2;
aB1 = @(pv) -fB1.beta ./ (pv - fB1.alpha).^2;

% En primer lugar, sin tener en cuenta los límites
%Dave = @(R, G, B) (fR1inv(R) + fG1inv(G) + fB1inv(B)) ./ 3;
% Y ahora teniéndolos en cuenta
isR = @(G) double(fG1inv(G)<10);
isB = @(G) double(fG1inv(G)>2);

Dave = @(R, G, B) (fR1inv(R).*isR(G) + fG1inv(G) + fB1inv(B).*isB(G)) ./ (isB(G) + isR(G) + 1);
wt = @(R,G,B) (aR1(R).*isR(G) + aG1(G) + aB1(B).*isB(G)).^2 ./ (isB(G) + isR(G) + 1) ./ ( (aR1(R).*isR(G)).^2 + (aG1(G)).^2 + (aB1(B).*isB(G)).^2 );

Dij = @(R,G,B) ( Dave(R,G,B) - wt(R,G,B).*(fR1inv(R).*aR1(R).*isR(G) + fG1inv(G).*aG1(G) + fB1inv(B).*aB1(B).*isB(G))./(aR1(R).*isR(G) + aG1(G) + aB1(B).*isB(G)) ) ./ (1 - wt(R,G,B));
figure(4);
plot(dosesGy, Dij(R,G,B), 'mo');
plot(dosesGy, Dave(R,G,B), 'cx');
resD = rssq((dosesGy'-Dij(R,G,B)), 1) / 10;
resDave = rssq((dosesGy'-Dave(R,G,B)), 1) / 10;
text(2,7,sprintf('res = %3.2f Gy',resD),'Color','magenta', 'FontSize', 14)
text(2,6,sprintf('res = %3.2f Gy',resDave),'Color','cyan', 'FontSize', 14)
legend({'Reference','Red channel', 'Green channel', 'Blue channel', 'Dij', 'Average D'}, 'Location', 'SouthEast');

%% Estudio sobre las imágenes completas (fórmula de la Derivada)
autoDoseValues = nan(Nsamples,1 );
autoDoseErrValues = nan(Nsamples, 1);
for i=1:Nsamples
    
    dosis1 = Dave(double(imageSubsets{i,filmIndex,1}(:,:,1)) ./ maxBits, ...
        double(imageSubsets{i,filmIndex,1}(:,:,2)) ./ maxBits, ...
        double(imageSubsets{i,filmIndex,1}(:,:,3)) ./ maxBits);
    dosis2 = Dave(double(imageSubsets{i,filmIndex,1}(:,:,1)) ./ maxBits, ...
        double(imageSubsets{i,filmIndex,1}(:,:,2)) ./ maxBits, ...
        double(imageSubsets{i,filmIndex,1}(:,:,3)) ./ maxBits);
    dosis3 = Dave(double(imageSubsets{i,filmIndex,1}(:,:,1)) ./ maxBits, ...
        double(imageSubsets{i,filmIndex,1}(:,:,2)) ./ maxBits, ...
        double(imageSubsets{i,filmIndex,1}(:,:,3)) ./ maxBits);    
    dosis = (dosis1+dosis2+dosis3)/3;
	figure(5);
    autoDoseValues(i) = mean2(dosis);
    autoDoseErrValues(i) = std2(dosis);
    imagesc(dosis)
    title(i);
    caxis([0 12]);
    colorbar
    fprintf('%3.3f +- %3.3f\n',mean2(dosis),std2(dosis));
    pause
end
figure
errorbar(dosesGy, autoDoseValues, autoDoseErrValues, 'b.');

%% Estudio sobre las imágenes completas (iterativa)
autoDoseValues = nan(Nsamples,1 );
autoDoseErrValues = nan(Nsamples, 1);
deltas = -0.2:0.002:0.2;

for i=1:Nsamples    
    I1a = imageSubsets{i,filmIndex,1};
    I1b = fliplr(I1a);
    I1c = flipud(I1a);
    I1d = fliplr(I1c);
    
    I2a = imageSubsets{i,filmIndex,2};
    I2b = fliplr(I2a);
    I2c = flipud(I2a);
    I2d = fliplr(I2c);
    
    I3a = imageSubsets{i,filmIndex,3};
    I3b = fliplr(I3a);
    I3c = flipud(I3a);
    I3d = fliplr(I3c);
    
    I = (double(I1a) + double(I1b) + double(I1c) + double(I1c) ...
      + double(I2a) + double(I2b) + double(I2c) + double(I2c) ...
      + double(I3a) + double(I3b) + double(I3c) + double(I3d)) / 12;
    pxmax = double(intmax(class(I1a)));
    
%     % Low-pass filtering code (NOT WORKING) --------------
%     
%     for RGBi=1:3
%         IR=I(:,:,RGBi); % convert the image to grey
%         A = fft2(double(Igrey)); % compute FFT of the grey image
%         Y=fftshift(A); % frequency scaling
%         [M N]=size(A); % image size
%         D=30; % filter size parameter
%         Lo(1:M,1:N)=0;
%         Lo(0.5*M-D:0.5*M+D,0.5*N-D:0.5*N+D)=1;
%         Lo = Lo(1:M, 1:N);
%         
%         % Filtered image=filter response*fft(original image)
%         J=Y.*Lo;
%         J1=ifftshift(J);
%         B1(:,:,RGBi)=real(ifft2(J1));
%     end
%     
%     % End low-pass filtering code ------------ 
    
    [dosis, varMat] = getDoseMicke(B1, CoefR, CoefG, CoefB, pixelsXcm, deltas, pxmax);
	figure(10+i);
    autoDoseValues(i) = mean2(dosis.data);
    autoDoseErrValues(i) = std2(dosis.data);
    subplot(1,2,1);
    dosis.plotSlice
    title(i);
    caxis([0 12]);
    colorbar
    subplot(1,2,2);
    varMat.plotSlice  
    caxis([min(deltas) max(deltas)]);
    colorbar
    fprintf('%3.3f +- %3.3f\n',autoDoseValues(i),autoDoseErrValues(i));
end
figure
errorbar(dosesGy, autoDoseValues, autoDoseErrValues, 'b.');
