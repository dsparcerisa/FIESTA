clear all
%%Cargar tablas de los polinomios
load('PolynomialEnergyDeposited.mat'); 
load('PolynomialFluence.mat');
%%Convertir las tablas en matrices
PFLU=[PolynomialFluence.E,PolynomialFluence.Afluence,PolynomialFluence.Bfluence,PolynomialFluence.Cfluence];
PEDEP=[PolynomialEnergyDeposited.E,PolynomialEnergyDeposited.AEdep,PolynomialEnergyDeposited.BEdep,PolynomialEnergyDeposited.CEdep];
E=5.3; %Energía a interpolar
Z=30;  %Plano
[Sigmaflu,SigmaEdep]=SigmaZE(PFLU,PEDEP,Z,E)