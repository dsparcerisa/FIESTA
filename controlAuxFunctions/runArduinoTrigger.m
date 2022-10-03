function arduinoStatus = runArduinoTrigger(puerto_serial, aux_npul, aux_width, aux_dc)

tic
%%Save the irradiation conditions in a log file 
% acqfolder=strcat('acq/',datestr(now,'yyyymmdd_HHMM'));
% mkdir(acqfolder);
% logfile=fopen(strcat(acqfolder,'/info.txt'),'w');
% log_npul=strcat('# pulses = ',aux_npul);
% log_width=strcat('Width (us) = ',aux_width);
% log_dc=strcat('Width Off (us) = ',aux_dc);
% 
% fprintf(logfile,'%s\n\r',log_npul);
% fprintf(logfile,'%s\n\r',log_width);
% fprintf(logfile,'%s\n\r',log_dc);


parameters = sprintf('%i,%i,%i,',aux_npul,aux_dc,aux_width);
display(parameters);

fprintf(puerto_serial, parameters);
pause(0.010)
nPulses=fscanf(puerto_serial,'%u')
widthOff_f=fscanf(puerto_serial,'%f')
widthOn_f=fscanf(puerto_serial,'%f')
widthOff_i=fscanf(puerto_serial,'%u')
widthOn_i=fscanf(puerto_serial,'%u')

  %  nb=1;
    while 1==1
        nb=puerto_serial.BytesAvailable;
        if nb>12
            arduinoStatus=fscanf(puerto_serial,'%s',nb);
            break;
        end
        pause(0.01)
    end
arduinoStatus
%log_widthOff=strcat('Width_Off (read from Arduino) (us) = ',widthOff);
%fprintf(logfile,'%s\n\r',log_widthOff);
%fclose(logfile);
toc
end