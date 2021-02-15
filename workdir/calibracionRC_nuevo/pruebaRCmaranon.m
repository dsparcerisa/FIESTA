clear all
filePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_10_13 Maranon/Maranon_1.tif';
I = imread(filePath);
I_irradiada = imcrop(I);
I_noirradiada = imcrop(I);
[pixCM, maxBits] = getImgMetaInfo(filePath);

I_R_irr = double(I_irradiada(:,:,1));
R_irr = mean(I_R_irr(:)) / maxBits
dR_irr = std(I_R_irr(:)) / maxBits

I_G_irr = double(I_irradiada(:,:,2));
G_irr = mean(I_G_irr(:)) / maxBits
dG_irr = std(I_G_irr(:)) / maxBits

I_B_irr = double(I_irradiada(:,:,3));
B_irr = mean(I_B_irr(:)) / maxBits
dB_irr = std(I_B_irr(:)) / maxBits

I_R_noirr = double(I_noirradiada(:,:,1));
R_noirr = mean(I_R_noirr(:)) / maxBits
dR_noirr = std(I_R_noirr(:)) / maxBits

I_G_noirr = double(I_noirradiada(:,:,2));
G_noirr = mean(I_G_noirr(:)) / maxBits
dG_noirr = std(I_G_noirr(:)) / maxBits

I_B_noirr = double(I_noirradiada(:,:,3));
B_noirr = mean(I_B_noirr(:)) / maxBits
dB_noirr = std(I_B_noirr(:)) / maxBits

%% Pruebas de dosis
load('fitCoef_HF_EBT3.mat')
applyRatFit = @(d, Coef) (Coef(2).*d + Coef(1)) ./ (d + Coef(3));
solveDose = @(PV, Coef) (Coef(1) - PV.*Coef(3)) ./ (PV - Coef(2));

DR = getChannelDoseT1(CoefR1, R_irr)
DG = getChannelDoseT1(CoefG1, G_irr)
DB = getChannelDoseT1(CoefB1, B_irr)


DRn = getChannelDoseT1(CoefR1, R_noirr)
DGn = getChannelDoseT1(CoefG1, G_noirr)
DBn = getChannelDoseT1(CoefB1, B_noirr)

DR2 = DR-DRn
DG2 = DG-DGn
DB2 = DB-DBn
D2 = (DR2 + DG2 + DB2) / 3

%% NetOD
D_R_NOD = getChannelDoseT3(CoefR3, R_irr, R_noirr)
D_G_NOD = getChannelDoseT3(CoefG3, G_irr, G_noirr)
D_B_NOD = getChannelDoseT3(CoefB3, B_irr, B_noirr)
D_NOD = (D_R_NOD + D_G_NOD + D_B_NOD) / 3

%% Micke aditivo
deltas = -0.2:0.001:0.2;
[dose, varMat] = getDoseMicke(double(I_irradiada), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
[dose_noIrr, varMat_noIrr] = getDoseMicke(double(I_noirradiada), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
D_Madd = mean2(dose.data) - mean2(dose_noIrr.data)

%% Micke multiplicativo
[dose2, varMat] = getDoseMicke_mult(double(I_irradiada), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
[dose2_noIrr, varMat_noIrr] = getDoseMicke_mult(double(I_noirradiada), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
D_Madd2 = mean2(dose2.data) - mean2(dose2_noIrr.data)


