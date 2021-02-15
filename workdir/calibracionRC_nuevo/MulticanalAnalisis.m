clc
clear all
close all

load CoefFit-fotones-EBT3unlnew.mat
[filename, filepath] = uigetfile('../EBT3unl/*.tif', 'Selecciona una imagen clinica'); 

Imagen = [filepath, filename]
I = imread(Imagen); 
ROI = imcrop(I);
pix_size_scan = 2.56/72; % 72 pix/inch, 1 inch = 2.56 cm
fscale=pix_size_scan/0.1; % resize pix_size from 0.0356cm to 0.1 cm
Image=imresize(ROI,fscale,'bicubic'); % resize image with a bicubic interpolation

dev = zeros(size(Image,1),size(Image,2));
dev3 = zeros(size(Image,1),size(Image,2));
DR = zeros(size(Image,1),size(Image,2));
DG = zeros(size(Image,1),size(Image,2));
DB = zeros(size(Image,1),size(Image,2));

figure(3); imagesc(Image); title('Clinical image');

ImageODR = log10(65535./double(Image(:,:,1)));
ImageODG = log10(65535./double(Image(:,:,2)));
ImageODB = log10(65535./double(Image(:,:,3)));

% Perturbed dose 
delta = 0.5:0.01:1.5;
dev_min = 1.e30*ones(size(Image,1),size(Image,2));

for k = 1:size(delta,2)
    for i = 1:size(Image,1)
        for j = 1:size(Image,2)
            DR(i,j) = (CoefR(3)*10^(-ImageODR(i,j)*delta(k))-CoefR(1))/(CoefR(2)-10^(-ImageODR(i,j)*delta(k)));
            DG(i,j) = (CoefG(3)*10^(-ImageODG(i,j)*delta(k))-CoefG(1))/(CoefG(2)-10^(-ImageODG(i,j)*delta(k)));
            DB(i,j) = (CoefB(3)*10^(-ImageODB(i,j)*delta(k))-CoefB(1))/(CoefB(2)-10^(-ImageODB(i,j)*delta(k)));
            
            dev(i,j) = (DR(i,j)-DG(i,j))^2+(DR(i,j)-DB(i,j))^2+(DB(i,j)-DG(i,j))^2;
            
            if dev(i,j) < dev_min(i,j)
                dev_min(i,j) = dev(i,j);
                delta0(i,j) = delta(k);
            end
        end
    end
    DR = zeros(size(Image,1),size(Image,2)); DG = DR; DB = DR;
end

DR = (CoefR(3)*10.^(-ImageODR.*delta0)-CoefR(1))./(CoefR(2)-10.^(-ImageODR.*delta0));
DG = (CoefG(3)*10.^(-ImageODG.*delta0)-CoefG(1))./(CoefG(2)-10.^(-ImageODG.*delta0));
DB = (CoefB(3)*10.^(-ImageODB.*delta0)-CoefB(1))./(CoefB(2)-10.^(-ImageODB.*delta0));

dev = (DR-DG).^2+(DR-DB).^2+(DB-DG).^2;
dev3 = (DR+DG+DB)./3;
dev2 = (DR+DG)./2;
DR_DB = DR./DB;

figure(4); imagesc(delta0); title('Perturbation (delta0)');
%figure(5); imagesc(dev_min); title('Deviation (dev_min)');
%figure(6); subplot(2,2,1);imagesc(dev3); title('Dose map (dev3)'); subplot(2,2,2); imagesc(DR); title('Dose red canal'); subplot(2,2,3); imagesc(DG); title('Dose green canal'); subplot(2,2,4); imagesc(DB); title('Dose blue canal'); 
%figure(6); subplot(2,2,1);imagesc(dev2); title('Dose map (dev2)'); subplot(2,2,2); imagesc(DR); title('Dose red canal'); subplot(2,2,3); imagesc(DG); title('Dose green canal');

%figure(7); subplot(2,1,1); imagesc(dev2-DR); title('Dose deviation red canal'); subplot(2,1,2); imagesc(dev2-DG); title('Dose deviation green canal');

figure(8); imagesc(dev3); title('Dose map (dev3)');

% Write converted dose in a text file

pix_size = 0.1; % pixel=0.1x0.1cm
fid=fopen([filename '.dat'],'w');

for i=1:size(dev3,1)
    for j=1:size(dev3,2)
       x = i*pix_size;
       y = j*pix_size;
       dosis = dev3(i,j);
       fprintf(fid,'%f %f %f\n',x,y,dosis);
    end
end

 fclose(fid);

% Reference doses

Dose2 = mean(mean(dev2,2))
Dose3 = mean(mean(dev3,2))
% Std = 100*(Dose - 4.0268)/4.0268; % in %
