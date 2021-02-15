function [dev3,delta0] = getDosePaula(Image,CoefR,CoefG,CoefB,pxmax)
ImageODR = log10(pxmax./double(Image(:,:,1)));
ImageODG = log10(pxmax./double(Image(:,:,2)));
ImageODB = log10(pxmax./double(Image(:,:,3)));

% Perturbed dose
delta = 0.95:0.002:1.05;
dev_min = 1.e30*ones(size(Image,1),size(Image,2));

for k = 1:numel(delta)
    DR = zeros(size(Image));
    DG = zeros(size(Image));
    DB = zeros(size(Image));
        
    for i = 1:size(Image,1)
        for j = 1:size(Image,2)
            Rij = 10^(-ImageODR(i,j)*delta(k));
            Gij = 10^(-ImageODR(i,j)*delta(k));
            Bij = 10^(-ImageODR(i,j)*delta(k));
            DR(i,j) = CoefR(3) + CoefR(2) / (Rij - CoefR(1));            
            DG(i,j) = CoefG(3) + CoefG(2) / (Gij - CoefG(1));            
            DB(i,j) = CoefB(3) + CoefB(2) / (Bij - CoefB(1));            
            
            dev(i,j) = (DR(i,j)-DG(i,j))^2+(DR(i,j)-DB(i,j))^2+(DB(i,j)-DG(i,j))^2;
            
            if dev(i,j) < dev_min(i,j)
                dev_min(i,j) = dev(i,j);
                delta0(i,j) = delta(k);
            end
        end
    end
end
DR = CoefR(3) + CoefR(2) ./ (10.^(-ImageODR.*delta0) - CoefR(1));            
DG = CoefG(3) + CoefG(2) ./ (10.^(-ImageODG.*delta0) - CoefG(1));  
DB = CoefB(3) + CoefB(2) ./ (10.^(-ImageODB.*delta0) - CoefB(1));  

dev3 = (DR+DG+DB)./3;
end

