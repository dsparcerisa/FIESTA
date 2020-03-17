function simRC = simulateRC(dose, Fcfg)
load(Fcfg.filmCalFile);
R = (CoefR(2)*dose + CoefR(1)) ./ (CoefR(3) + dose);
G = (CoefG(2)*dose + CoefG(1)) ./ (CoefG(3) + dose);
B = (CoefB(2)*dose + CoefB(1)) ./ (CoefB(3) + dose);
Ri = uint8(round(255*R));
Gi = uint8(round(255*G));
Bi = uint8(round(255*B));
simRC = uint8(zeros(size(Ri,1), size(Ri,2), 3));
simRC(:,:,1) = Ri;
simRC(:,:,2) = Gi;
simRC(:,:,3) = Bi;
end

