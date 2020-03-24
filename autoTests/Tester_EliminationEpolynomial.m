clear all

load('PolynomialFluence.mat')
load('PolynomialEnergyDeposited.mat')

%Set energy and range to test the interpolation polynomial error per plane

E=9;
Z=1:20;

Errorflu=nan(1,20);
ErrorEdep=nan(1,20);
for i=1:length(Z);
    [Sigmaflu,SigmaEdep]=EliminationEpolynomial(PolynomialFluence,PolynomialEnergyDeposited,Z(i),E);
    Errorflu(i)=Sigmaflu(3);
    ErrorEdep(i)=SigmaEdep(3);
end
varNames={'Z','ErrorFluence','ErrorEnergyDep'};
tblErrorPerPlane=table(Z',Errorflu',ErrorEdep','VariableNames',varNames)