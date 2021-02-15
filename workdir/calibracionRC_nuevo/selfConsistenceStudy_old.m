%% Estudio de la solución de rat12
%syms pv p1 D p2 q1 q2
%rat12 = pv == (p1*D + p2)/(D^2 + q1*D + q2)
%S = solve(rat12, D)

fR1inv = @(pv) fR1.gamma + fR1.beta./(pv - fR1.alpha);
fG1inv = @(pv) fG1.gamma + fG1.beta./(pv - fG1.alpha);
fB1inv = @(pv) fB1.gamma + fB1.beta./(pv - fB1.alpha);

fB1binv_aux = @(pv,p1,p2,q1,q2) (p1 + (p1.^2 - 2.*p1.*pv.*q1 + pv.^2.*q1.^2 - 4.*q2.*pv.^2 + 4.*p2.*pv).^(1/2) - pv.*q1)./(2.*pv)
real([dosesGy' fR1inv(R') fG1inv(G') fB1inv(B')])

%% Estudio de autoconsistencia para todos:
figure(4)
dosesR = fR1inv(R);
ddosesR = max(abs(fR1inv(R+dR)-dosesR), abs(fR1inv(R-dR)-dosesR));
dosesG = fG1inv(G);
ddosesG = max(abs(fG1inv(G+dG)-dosesG), abs(fG1inv(G-dG)-dosesG));
dosesB = fB1inv(B);
ddosesB = max(abs(fB1inv(B+dB)-dosesB), abs(fB1inv(B-dB)-dosesB));

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
resR = rssq((dosesGy-fR1inv(R))) / Nsamples;
resG = rssq((dosesGy-fG1inv(G))) / Nsamples;
resB = rssq((dosesGy-fB1inv(B))) / Nsamples;
axis([0 12 0 12]);

text(2,10,sprintf('res = %3.2f Gy',resR),'Color','red', 'FontSize', 14)
text(2,9,sprintf('res = %3.2f Gy',resG),'Color','green', 'FontSize', 14)
text(2,8,sprintf('res = %3.2f Gy',resB),'Color','blue', 'FontSize', 14)
set(gca, 'FontSize', 14)

%% Estudio con el Mayer para Tipo I
aR1 = @(pv) -fR1.beta ./ (pv - fR1.alpha).^2;
aG1 = @(pv) -fG1.beta ./ (pv - fG1.alpha).^2;
aB1 = @(pv) -fB1.beta ./ (pv - fB1.alpha).^2;

% En primer lugar, sin tener en cuenta los límites
%Dave = @(R, G, B) (fR1inv(R) + fG1inv(G) + fB1inv(B)) ./ 3;
% Y ahora teniéndolos en cuenta
%isR = @(G) double(fG1inv(G)<10);
%isB = @(G) double(fG1inv(G)>2);
isR = @(G) 1.0;
isB = @(G) 1.0;

Dave = @(R, G, B) (fR1inv(R).*isR(G) + fG1inv(G) + fB1inv(B).*isB(G)) ./ (isB(G) + isR(G) + 1);
wt = @(R,G,B) (aR1(R).*isR(G) + aG1(G) + aB1(B).*isB(G)).^2 ./ (isB(G) + isR(G) + 1) ./ ( (aR1(R).*isR(G)).^2 + (aG1(G)).^2 + (aB1(B).*isB(G)).^2 );

Dij = @(R,G,B) ( Dave(R,G,B) - wt(R,G,B).*(fR1inv(R).*aR1(R).*isR(G) + fG1inv(G).*aG1(G) + fB1inv(B).*aB1(B).*isB(G))./(aR1(R).*isR(G) + aG1(G) + aB1(B).*isB(G)) ) ./ (1 - wt(R,G,B));
figure(4);
plot(dosesGy, dosesGy, 'k-'); hold on
plot(dosesGy, fR1inv(R), 'ro');
plot(dosesGy, fG1inv(G), 'go');
plot(dosesGy, fB1inv(B), 'bo');
plot(dosesGy, Dij(R,G,B), 'mo'); 
plot(dosesGy, Dave(R,G,B), 'cx');
resD = rssq(dosesGy-Dij(R,G,B)) / Nsamples;
resDave = rssq(dosesGy-Dave(R,G,B)) / Nsamples;

text(1,11,sprintf('res = %3.2f Gy',resR),'Color','red', 'FontSize', 14)
text(1,10,sprintf('res = %3.2f Gy',resG),'Color','green', 'FontSize', 14)
text(1,9,sprintf('res = %3.2f Gy',resB),'Color','blue', 'FontSize', 14)
text(1,8,sprintf('res = %3.2f Gy',resD),'Color','magenta', 'FontSize', 14)
text(1,7,sprintf('res = %3.2f Gy',resDave),'Color','cyan', 'FontSize', 14)
ylim([0 round(max(dosesGy))+2]);
grid on;
legend({'Reference','Red channel', 'Green channel', 'Blue channel', 'Dij (Mayer formula)', 'Average D'}, 'Location', 'SouthEast');

%% Estudio sobre las imágenes completas (Dave, promedio de las tres)
autoDoseValues = nan(Nsamples,1 );
autoDoseErrValues = nan(Nsamples, 1);
for j=1:Nsamples
    
    dosis1 = Dave(double(imageSubsets{filmIndex,j,1}(:,:,1)) ./ maxBits, ...
        double(imageSubsets{filmIndex,j,1}(:,:,2)) ./ maxBits, ...
        double(imageSubsets{filmIndex,j,1}(:,:,3)) ./ maxBits);
    dosis2 = Dave(double(imageSubsets{filmIndex,j,2}(:,:,1)) ./ maxBits, ...
        double(imageSubsets{filmIndex,j,2}(:,:,2)) ./ maxBits, ...
        double(imageSubsets{filmIndex,j,2}(:,:,3)) ./ maxBits);
    dosis3 = Dave(double(imageSubsets{filmIndex,j,3}(:,:,1)) ./ maxBits, ...
        double(imageSubsets{filmIndex,j,3}(:,:,2)) ./ maxBits, ...
        double(imageSubsets{filmIndex,j,3}(:,:,3)) ./ maxBits);    
    dosis = (dosis1+dosis2+dosis3)/3;
	figure(5);
    autoDoseValues(j) = mean2(dosis);
    autoDoseErrValues(j) = std2(dosis);
    %imagesc(dosis)
    %title(i);
    %caxis([0 12]);
    %colorbar
    fprintf('%3.3f +- %3.3f\n',mean2(dosis),std2(dosis));
    %pause
end
figure
errorbar(dosesGy, autoDoseValues, autoDoseErrValues, 'b.');
grid on
xlim([0 round(max(dosesGy))]);
ylim([0 round(max(dosesGy))+2]);

%% Estudio sobre las imágenes completas (iterativa)
autoDoseValues = nan(Nsamples,1);
autoDoseErrValues = nan(Nsamples,1);
deltas = -0.2:0.002:0.2;

for j=1:Nsamples    
    I1a = imageSubsets{filmIndex,j,1};
    I1b = fliplr(I1a);
    I1c = flipud(I1a);
    I1d = fliplr(I1c);
    
    I2a = imageSubsets{filmIndex,j,2};
    I2b = fliplr(I2a);
    I2c = flipud(I2a);
    I2d = fliplr(I2c);
    
    I3a = imageSubsets{filmIndex,j,3};
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
    
    [dosis, varMat] = getDoseMicke(I, CoefR, CoefG, CoefB, pixelsXcm, deltas, pxmax);
	figure(10+j);
    autoDoseValues(j) = mean2(dosis.data);
    autoDoseErrValues(j) = std2(dosis.data);
    subplot(1,2,1);
    dosis.plotSlice
    title(j);
    caxis([0 round(max(dosesGy))]);
    colorbar
    subplot(1,2,2);
    varMat.plotSlice  
    caxis([min(deltas) max(deltas)]);
    colorbar
    fprintf('%3.3f +- %3.3f\n',autoDoseValues(j),autoDoseErrValues(j));
end
figure
errorbar(dosesGy, autoDoseValues, autoDoseErrValues, 'b.');
grid on