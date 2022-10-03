% Script to start instruments

% 1. Connect

if  ~exist('COMtrigger','var') || strcmp(COMtrigger.Status,'closed')
    COMtrigger = serial(Fcfg.triggerPort,'BaudRate',9600);
    
    fopen(COMtrigger); disp('Trigger connected');       
end