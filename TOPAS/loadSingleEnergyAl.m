function [Sflu, SEdep, DoseZ] = loadSingleEnergy(depthString)
% Creates polinomial fit for sigma based on fluence and deposited energy
% in files within folder

%% Prepare file names
EDepFileName = 'EnergyDep.csv';
FluenceFileName = 'Fluence.csv';
%% Set voxel dimensions
% # Results for scorer Fluence
% # R in 100 bins of 0.01 cm
% # Phi in 1 bin  of 360 deg
% # Z in 2000 bins of 0.02 cm
% # Fluence ( /mm2 ) : Sum
% # Edep ( MeV ) : Sum

NR = 300;
NZ = 2000;
dR = 0.01;
dZ = 0.02;

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [9, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["A", "B", "C", "D"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

%% Read files for desired energy
tblFlu = readtable(fullfile('AlFoil', depthString, FluenceFileName), opts);
tblEdep = readtable(fullfile('AlFoil', depthString, EDepFileName), opts);

% Create fluence and Edep matrices
flu = reshape(tblFlu.D, [NZ NR]);
Edep = reshape(tblEdep.D, [NZ NR]);

%% Calculate and plot depth dose distribution
rho = 1;
RMaxValues = dR * (1:NR);
RMinValues = RMaxValues - dR;
RValues = RMaxValues - dR/2;
ZMaxValues = dZ * (1:NZ);
ZMinValues = ZMaxValues - dZ;
ZValues = ZMaxValues - dZ/2;

A = pi*(RMaxValues.^2-RMinValues.^2);
Atot = sum(A);
Dose = Edep.* ((rho*dZ*A).^(-1));
DoseZ = sum(Dose.*A,2); %/Atot; % Total Dose per plane
DoseZ = DoseZ / sum(DoseZ);
EdepZ = sum(Edep,2);
EdepZ = EdepZ * DoseZ(1) / EdepZ(1);

figure;
subplot(2,1,1);
plot(ZValues, DoseZ, 'b-');
hold on
plot(ZValues, EdepZ, 'r:');
legend('Geometrically calculated dose', 'Transversal sum of Edep');
xlabel('Depth (cm)');
ylabel('Integrated dose (a.u.)');
title(depthString);

subplot(2,1,2);
imagesc(ZValues,RValues,Edep');
xlabel('Depth (cm)');
ylabel('Radius (cm)');
title('Edep');
axis xy
set(gca, 'colorscale', 'log');

% --> A partir de 5 MeV falla porque "no cabe" todo el haz dentro del scorer
% --> Habrá que recalcular las tablas con un scorer más ancho.

%% Find maximum value ("range");
rangeIndex = find(EdepZ==0,1);
if isempty(rangeIndex)
    rangeIndex = NZ;
else
    rangeIndex = rangeIndex - 1;
end

%% Sigma for the Fluence
Sflu = nan(1, NZ);
for i=1:rangeIndex
    try
        Fvalues = flu(i,:);
        Fvalues = Fvalues / sum(Fvalues);
        F1 = fit(RValues', Fvalues', 'gauss1');
        Sflu(i) = F1.c1 / sqrt(2);
    catch
        % minValidIndex = i;
    end
end
[~, maxPos] = max(Sflu);
maxValidIndex = min(maxPos, rangeIndex);

[~, minPos] = min(Sflu(1:round(maxValidIndex/2)));
minValidIndex = minPos;

%Sflu(1:minValidIndex) = nan;
%Sflu(maxValidIndex:NZ) = nan;

%% Sigma for the Edep
SEdep = nan(1, NZ);
for i=1:maxValidIndex
    try
        Fvalues = Edep(i,:);
        Fvalues = Fvalues / sum(Fvalues);
        F1 = fit(RValues', Fvalues', 'gauss1');
        SEdep(i) = F1.c1 / sqrt(2);
    catch
        % minValidIndex = i;
    end
end

%SEdep(1:minValidIndex) = nan;
%SEdep(maxValidIndex:NZ) = nan;


end

