% Script to start instruments

% 1. Connect

if  ~exist('COMshutter') || strcmp(COMshutter.Status,'closed')
    COMshutter = serial(Fcfg.arduinoPort,'BaudRate',9600);
    
    fopen(COMshutter); disp('Shutter connected');       
end