function [dose, varMat] = getDoseFromRC(I, CoefR, CoefG, CoefB, pixelsXCM, deltas)

pxmax = double(intmax(class(I))); % maximum pixel value 2^16

IODR = log10(pxmax./double(I(:,:,1)));
IODG = log10(pxmax./double(I(:,:,2)));
IODB = log10(pxmax./double(I(:,:,3)));

if ~exist('deltas')
    delta = 0.8:0.002:1.2;
else
    delta = deltas;
end

% Perturbed dose 
dev_min = 1.e30*ones(size(I,1),size(I,2));
delta0 = nan(size(I,1),size(I,2));

for k = 1:numel(delta)
    for i = 1:size(I,1)
        for j = 1:size(I,2)
            DR = (CoefR(3)*10^(-IODR(i,j)*delta(k))-CoefR(1))/(CoefR(2)-10^(-IODR(i,j)*delta(k)));
            DG = (CoefG(3)*10^(-IODG(i,j)*delta(k))-CoefG(1))/(CoefG(2)-10^(-IODG(i,j)*delta(k)));
            DB = (CoefB(3)*10^(-IODB(i,j)*delta(k))-CoefB(1))/(CoefB(2)-10^(-IODB(i,j)*delta(k)));
            
            dev = (DR-DG)^2+(DR-DB)^2+(DB-DG)^2;
            
            if dev < dev_min(i,j)
                dev_min(i,j) = dev;
                delta0(i,j) = delta(k);
            end
        end
    end
end

DR = (CoefR(3)*10.^(-IODR.*delta0)-CoefR(1))./(CoefR(2)-10.^(-IODR.*delta0));
DG = (CoefG(3)*10.^(-IODG.*delta0)-CoefG(1))./(CoefG(2)-10.^(-IODG.*delta0));
DB = (CoefB(3)*10.^(-IODB.*delta0)-CoefB(1))./(CoefB(2)-10.^(-IODB.*delta0));

doseRGB = (DR+DG+DB)./3;

%figure(1); imagesc(doseRGB); title('Dose map (RGB)'); caxis([0 10]);
%figure(2); imagesc(delta0); title('Delta0'); caxis([0.8 1.2]);

dxy = 1/pixelsXCM;
NX = size(doseRGB, 1);
NY = size(doseRGB, 2);
sizeX = dxy*(NX);
sizeY = dxy*(NY);
dose = createEmptyCG2D(dxy, sizeX, sizeY);
dose.data = doseRGB;

varMat = createEmptyCG2D(dxy, sizeX, sizeY);
varMat.data = delta0;

end

