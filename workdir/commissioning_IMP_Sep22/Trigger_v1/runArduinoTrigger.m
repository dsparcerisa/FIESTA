function status = runArduinoTrigger(puerto_serial, aux_npul, aux_width, aux_dc)


%%Save the irradiation conditions in a log file 
acqfolder=strcat('acq/',datestr(now,'yyyymmdd_HHMM'));
mkdir(acqfolder);
logfile=fopen(strcat(acqfolder,'/info.txt'),'w');
log_npul=strcat('# pulses = ',aux_npul);
log_width=strcat('Width (us) = ',aux_width);
log_dc=strcat('Width Off (us) = ',aux_dc);

fprintf(logfile,'%s\n\r',log_npul);
fprintf(logfile,'%s\n\r',log_width);
fprintf(logfile,'%s\n\r',log_dc);


parameters = sprintf('%i,%i,%i,',aux_npul,aux_dc,aux_width);
display(parameters);

tic
fprintf(puerto_serial, parameters);
widthOff=fscanf(puerto_serial,'%s');
status=fscanf(puerto_serial,'%s');
toc
log_widthOff=strcat('Width_Off (read from Arduino) (us) = ',widthOff);
fprintf(logfile,'%s\n\r',log_widthOff);
fclose(logfile);
end