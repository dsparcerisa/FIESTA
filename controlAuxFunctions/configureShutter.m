function configureShutter(COM,A,value)
% set Shutter ON time (seconds) for normal irradiation
    status='';
    nb=1;
    while nb>0
        nb=COM.BytesAvailable;
        if nb>0
            status=fscanf(COM,'%c',nb);
        end
    end

    fprintf(1,'Configurating Shutter %f seconds\n',value)
    fprintf(COM,strcat(A,num2str(value)));
    
    status='';
    while numel(strfind(status,'D'))==0
        nb=COM.BytesAvailable;
        if nb>0
            status=fscanf(COM,'%c',nb);
        end
    end    
    
    disp('done configurating shutter');
end