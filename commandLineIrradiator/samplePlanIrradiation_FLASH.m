% Add path to consider 
stageControlPath = '../stageControlGUI';
addpath(stageControlPath);

% Start instruments
startInstruments

% Read plan
plan = readPlan('plan.txt');

%% Data necessary for irradiation
t_ms = 0.65; % shot equivalent time in ms

Ifactor = 0.81;
I_FC1_nA = 1.90; 

Q_shot_pC = I_FC1_nA * Ifactor * t_ms;

plan.n_shots = round(plan.Q ./ Q_shot_pC);

%% Center stage
autoCenter(s2);   
[~, ~, ~, Zpos] = monitorStatus(s2);

%% Align in Z position
% align manually and then set here
Zpos = 0;

%% Irradiate plan
tic
for i=1:numel(plan.X)
    absPos = [plan.X(i) plan.Y(i) Zpos];
    finished = stageControl_moveToAbsPos(s2, absPos);
    Shutter(s1,'f',plan.n_shots(i));
    fprintf('Irradiated spot %i\n', i);
end
toc