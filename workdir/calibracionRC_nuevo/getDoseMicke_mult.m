function [dose, varMat, realDR, realDG, realDB] = getDoseMicke_mult(I, CoefR, CoefG, CoefB, pixelsXCM, deltas, pxmax)


%% Definir funciones
DR = @(pv) CoefR(3) + CoefR(2)./(pv - CoefR(1));
DG = @(pv) CoefG(3) + CoefG(2)./(pv - CoefG(1));
DB = @(pv) CoefB(3) + CoefB(2)./(pv - CoefB(1));
aR = @(pv) -CoefR(2) ./ (pv - CoefR(1)).^2 .* log(10) .* log10(pv);
aG = @(pv) -CoefG(2) ./ (pv - CoefG(1)).^2 .* log(10) .* log10(pv);
aB = @(pv) -CoefB(2) ./ (pv - CoefB(1)).^2 .* log(10) .* log10(pv);

correctedDR = @(pv,delta) DR(pv) + aR(pv).*delta;
correctedDG = @(pv,delta) DG(pv) + aG(pv).*delta;
correctedDB = @(pv,delta) DB(pv) + aB(pv).*delta;

pvR = double(I(:,:,1)) / pxmax;
pvG = double(I(:,:,2)) / pxmax;
pvB = double(I(:,:,3)) / pxmax;

if ~exist('deltas')
    delta = -0.2:0.001:0.2;
else
    delta = deltas;
end

% Perturbed dose 
dev_min = 1.e30*ones(size(I,1),size(I,2));
delta0 = nan(size(I,1),size(I,2));

% Primero sin tener en cuenta los límites
for k = 1:numel(delta)
    dev = (correctedDR(pvR,delta(k))-correctedDG(pvG,delta(k))).^2 + ...
        (correctedDR(pvR,delta(k))-correctedDB(pvB,delta(k))).^2 + ...
        (correctedDB(pvB,delta(k))-correctedDG(pvG,delta(k))).^2;
    maskMin = dev<dev_min;
    dev_min(maskMin) = dev(maskMin);
    delta0(maskMin) = delta(k);
end

realDR = correctedDR(pvR,delta0);
realDG = correctedDG(pvG,delta0);
realDB = correctedDB(pvB,delta0);

% Y luego teniéndolos en cuenta
% isR = double(DG(pvG)>0);
% isB = double(DG(pvG)>0);
% 
% for k = 1:numel(delta)
%     dev = isR.*(correctedDR(pvR,delta(k))-correctedDG(pvG,delta(k))).^2 + ...
%         isR.*isB.*(correctedDR(pvR,delta(k))-correctedDB(pvB,delta(k))).^2 + ...
%         isB.*(correctedDB(pvB,delta(k))-correctedDG(pvG,delta(k))).^2;
%     maskMin = dev<dev_min;
%     dev_min(maskMin) = dev(maskMin);
%     delta0(maskMin) = delta(k);
% end
% 
% doseRGB = (isR.*realDR+realDG+isB.*realDB)./ (1 + isR + isB);
doseRGB = (realDR + realDG + realDB) / 3;

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


