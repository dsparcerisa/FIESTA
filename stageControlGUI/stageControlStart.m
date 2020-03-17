function stageControlStart(COM)
if ~exist('COM')
    instrreset
    %% Connect
    % if ~exist('s1')
    %     s1 = serial(arduinoport,'BaudRate',9600);
    % end
    
    if  ~exist('s2') ||  strcmp(s2.Status,'closed')
        instrreset
        
        % s1 = serial(arduinoport,'BaudRate',9600);
        s2 = serial(Fcfg.velmexPort,'BaudRate',9600);
        
        % fopen(s1); disp('Shutter connected');
        fopen(s2); disp('Stages connected');
        
        fprintf(s2,'V');
        pause(.1); % wait for 100ms
        
        % see if the controller connected properly
        readStatus(s2);
        
        % Initialize stage
        disp('Initializing at speed 6000 and acceleration 50')
        fprintf(s2,'F,C,setM1M3,S1M6000,A1M50,setL1M1,R');
        fprintf(s2,'F,C,setM2M3,S2M6000,A2M50,setL2M1,R');
        fprintf(s2,'F,C,setM3M3,S3M6000,A3M50,setL3M1,R');
        
        readStatus(s2);
        
        set(s2,'Timeout',30);
        
    end
else
    s2 = COM;
end
%% Open instance of stageControl
stageControl_exported.getInstance(s2);
end