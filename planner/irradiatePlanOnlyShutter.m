function irradiatePlanOnlyShutter(Fcfg, COMshutter, numCycles, t_open, t_wait)
% void irradiatePlan(COMstage, COMshutter, plan, vector2startPoint, stageLimits, airDepthAtPos0, app)
% Launches plan irradiation
if (exist('app'))
    app.logLine('Starting plan irradiation...');
else
    disp('Starting plan irradiation...');
end

%% Create log file
logFileName = datestr(now,'irrLog_yyyy_mm_dd_HH_MM_SS.log');
logFile = fopen(fullfile(Fcfg.logPath,logFileName),'w');

msg = sprintf('%i cycles of %3.3fs open and %3.3fs waiting', numCycles, t_open, t_wait);  
fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);

fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') 'Starting plan irradiation...' '\n']);
comment = input('Enter plan description text for log: ','s');
fprintf(logFile, '\tPlan comment: %s\n', comment);

%% Input intensity at FC1
I_nA = input('Input currrent intensity at FC1 (nA): ');
fprintf(logFile, '\tCurrent intensity: %3.3f nA \n\n', I_nA);

    
msg = sprintf('Configure shutter to open for for %3.3fs', t_open);  
fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
configureShutter(COMshutter, 't', t_open)
pause(0.1);
msg = sprintf('Done shutter configuration.\n');  
fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
    
%% Irradiate plan
for i=1:numCycles
    msg = sprintf('Starting cycle %i of %i', i, numCycles);       
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);

    msg = sprintf('Opening shutter for %3.3fs', t_open);  
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
    shutter(COMshutter,'n',1); 
    msg = sprintf('Shutter closed.\n');       
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
               
    fprintf('Irradiated cycle %i / %i\n', i, numCycles);
    pause(t_wait);
end


%% Close log file and end
msg = sprintf('Plan irradiation finished.');       
fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
fclose(logFile);

if (exist('app'))
    app.logLine('Plan irradiation finished.');
else
    disp('Plan irradiation finished.');
end

end

