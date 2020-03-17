% Barrer en XYZ

initPos = [0 0 0];

destPos = [3 3 3; 1 1 1; -1 -1 -1; -3 -3 -3; -5 -5 -5];
steps = 2:2:10;
N = 5;
timeXminus = nan(N,1);
timeXplus = nan(N,1);
for i = 1:size(destPos, 1)
    finished = stageControl_moveToAbsPos(s2, initPos);  
    tic
    finished = stageControl_moveToAbsPos(s2, destPos(i, :));
    timeXminus(i) = toc
    tic
    finished = stageControl_moveToAbsPos(s2, initPos);  
    timeXplus(i) = toc    
    toc
end

%% Plottear
plot(steps, timeXminus, 'bo'); hold on;
plot(steps, timeXplus, 'ro')
Fminus = fit(3*steps', timeXminus, 'poly1')
Fplus = fit(3*steps', timeXplus, 'poly1')
% Z Time = 1.89 + 0.66*d(cm);
% Y Time = 1.24 + 0.66*d(cm);
% X Time = 1.17 + 0.66*d(cm);
% XY Time = 2.03 + 0.66*d(cm);
% XZ Time = 1.13 + 0.66*d(cm);
% YZ Time = 1.11 + 0.66*d(cm);
% XYZ Time = 1.83 + 0.66*d(cm);

%% Realizar una prueba de un 3x3 en XY (sin shutter)
delta = 0.2; % cm
initPos = [-delta -delta 0];
deltaTimes = nan(8,1);
finished = stageControl_moveToAbsPos(s2, initPos);  

initTime = tic;

tic
initPos = [-delta 0 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(1) = toc;

tic
initPos = [-delta delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(2) = toc;

tic
initPos = [0 delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(3) = toc;

tic
initPos = [0 0 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(4) = toc;

tic
initPos = [0 -delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(5) = toc;

tic
initPos = [delta -delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(6) = toc;

tic
initPos = [delta 0 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(7) = toc;

tic
initPos = [delta +delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(8) = toc;

totaltime = toc(initTime)

% With delta = 0.1 da 9.93 (suma) o totalTime (da lo mismo) 
% --> Lag time 1.18
% With delta = 0.2 esperamos 10.5 (suma), da 9.93

%% Realizar una prueba de un 3x3 en XY (CON shutter)
delta = 0.2; % cm
timeIrr = 2;
Nshots = 1;
initPos = [-delta -delta 0];
deltaTimes = nan(8,1);
finished = stageControl_moveToAbsPos(s2, initPos);  
Configure_shutter(s1, 't', timeIrr);
Shutter(s1,'n',Nshots);
initTime = tic;

tic
initPos = [-delta 0 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(1) = toc;

tic
initPos = [-delta delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(2) = toc;

tic
initPos = [0 delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(3) = toc;

tic
initPos = [0 0 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(4) = toc;

tic
initPos = [0 -delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(5) = toc;

tic
initPos = [delta -delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(6) = toc;

tic
initPos = [delta 0 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(7) = toc;

tic
initPos = [delta +delta 0];
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(8) = toc;

totaltime = toc(initTime)

% With delta = 0.1 da 9.93 (suma) o totalTime (da lo mismo) 
% --> Lag time 1.18
% With delta = 0.2 esperamos 10.5 (suma)
% TIME = 9.75 + 2.67 * Nshots/spot [para Flash], sale 0.3s por shot
% TIME = 10.1 + 9 * tirr [para Normal], lo esperable.


%% Realizar una prueba de un 2x3x2 en XY (sin shutter)
delta = 0.2; % cm
initPos = [-delta/2 -delta*cosd(30) 0]; %1
deltaTimes = nan(6,1);
finished = stageControl_moveToAbsPos(s2, initPos);  
initTime = tic;

tic
initPos = [+delta/2 -delta*cosd(30) 0]; %2
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(1) = toc;

tic
initPos = [+delta/2 +delta*cosd(30) 0]; %3
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(2) = toc;
tic
initPos = [-delta/2 +delta*cosd(30) 0]; %4
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(3) = toc;

tic
initPos = [-delta 0 0]; %5
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(4) = toc;

tic
initPos = [0 0 0]; %6
finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(5) = toc;

tic
initPos = [delta 0 0]; %7

finished = stageControl_moveToAbsPos(s2, initPos); 
deltaTimes(6) = toc;

totaltime = toc(initTime)

% TOTAL TIME = 4.9 s (tras cambiar el pause en readStatus a 0.04 s)

%% Prueba de un 2x3x2 en XY (con shutter)
delta = 0.2; % cm
Nshots = 1;
irrT = 2;
initPos = [-delta/2 -delta*cosd(30) 0]; %1
deltaTimes = nan(6,1);
finished = stageControl_moveToAbsPos(s2, initPos);  
initTime = tic;

Configure_shutter(s1,'t',irrT);
pause(0.05);
Shutter(s1,'n',Nshots);

tic
initPos = [+delta/2 -delta*cosd(30) 0]; %2
finished = stageControl_moveToAbsPos(s2, initPos); 
Configure_shutter(s1,'t',irrT);
pause(0.05);
Shutter(s1,'n',Nshots);
deltaTimes(1) = toc;

tic
initPos = [+delta/2 +delta*cosd(30) 0]; %3
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(2) = toc;
tic
initPos = [-delta/2 +delta*cosd(30) 0]; %4
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(3) = toc;

tic
initPos = [-delta 0 0]; %5
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(4) = toc;

tic
initPos = [0 0 0]; %6
finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(5) = toc;

tic
initPos = [delta 0 0]; %7

finished = stageControl_moveToAbsPos(s2, initPos); 
Shutter(s1,'n',Nshots);
deltaTimes(6) = toc;

totaltime = toc(initTime)