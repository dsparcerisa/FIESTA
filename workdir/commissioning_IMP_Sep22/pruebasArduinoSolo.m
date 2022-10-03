%% Start trigger
clear COMtrigger
instrreset
FIESTAConfig_win_sep22
startTrigger
%% Test
nShots = 3
wIn = 10;
wOut = 1000000;
msg = runArduinoTrigger(COMtrigger,  nShots, wIn, wOut);      

%% Cerrar
fclose(COMtrigger)