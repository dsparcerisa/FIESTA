% Script to start Stage

if ~exist('varSpeed')
    varSpeed = 6000;
    warning('Variable varSpeed not specified. Using default max value of 6000');
end

% 1. Connect

    COMstage = serial(Fcfg.velmexPort,'BaudRate',9600);
    
    fopen(COMstage); disp('Stages connected');
    
    fprintf(COMstage,'V');
    pause(.1); % wait for 100ms
    
    % see if the controller connected properly
    readStatus(COMstage);
    
    % Initialize stage
    fprintf('Initializing at speed %i and acceleration 50', varSpeed)
    fprintf(COMstage,'F,C,setM1M3,S1M6000,A1M50,setL1M1,R');
    fprintf(COMstage,'F,C,setM2M3,S2M6000,A2M50,setL2M1,R');
    fprintf(COMstage,'F,C,setM3M3,S3M6000,A3M50,setL3M1,R');
    
    readStatus(COMstage);
    
    set(COMstage,'Timeout',30);
    