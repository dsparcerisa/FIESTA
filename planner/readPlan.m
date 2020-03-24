function plan = readPlan(planFile)

%% Read name
ID = fopen(planFile);
tline = fgetl(ID);
tline = strtrim(tline(2:end));
fprintf('Reading plan with name %s\n', tline);

%% Read plan type
planType = fgetl(ID);
fclose(ID);

if strcmp(planType, 'FLASH')
    disp('FLASH type plan');
elseif strcmp(planType, 'CONV')
    disp('CONV type plan');   
else
    error('Unrecognized plan type (%s)', planType);
end

%% Read header

opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [4, 4];
opts.Delimiter = ",";

% Specify column names and types
if strcmp(planType, 'FLASH')
    opts.VariableNames = ["EMeV", "InA", "tRendija", "numSpots"];
    opts.VariableTypes = ["double", "double", "double", "uint16"];    
elseif strcmp(planType, 'CONV')
    opts.VariableNames = ["EMeV", "InA", "codFiltro", "numSpots"];
    opts.VariableTypes = ["double", "double", "string", "uint16"];
    opts = setvaropts(opts, "codFiltro", "WhitespaceRule", "trim");
    opts = setvaropts(opts, "codFiltro", "EmptyFieldRule", "auto");
else
    error('Unrecognized plan type (%s)', planType);
end

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "numSpots", "TrimNonNumeric", true);
opts = setvaropts(opts, "numSpots", "ThousandsSeparator", ",");

% Import header
header = readtable(planFile, opts);
clear opts

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [6, Inf];
opts.Delimiter = ",";

% Specify column names and types

if strcmp(planType, 'FLASH')
    opts.VariableNames = ["Xcm", "Ycm", "Zcm", "Nshots"];
    opts.VariableTypes = ["double", "double", "double", "uint16"];    
elseif strcmp(planType, 'CONV')
    opts.VariableNames = ["Xcm", "Ycm", "Zcm", "Ts"];
    opts.VariableTypes = ["double", "double", "double", "double"];
else
    error('Unrecognized plan type (%s)', planType);
end


% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
dataTable = readtable(planFile, opts);

if strcmp(planType, 'FLASH')
    plan.tRendija = header.tRendija;
    plan.Nshots = dataTable.Nshots;        
elseif strcmp(planType, 'CONV')
    plan.codFiltro = header.codFiltro;
    plan.t_s = dataTable.Ts;    
else
    error('Unrecognized plan type (%s)', planType);
end

plan.E = header.EMeV;
plan.I = header.InA;
plan.numSpots = header.numSpots;
plan.X = dataTable.Xcm;
plan.Y = dataTable.Ycm;
plan.Z = dataTable.Zcm;
plan.name = tline;
plan.mode = planType;

end

