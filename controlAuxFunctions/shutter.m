function shutter(COM,A,i)

for n=1:i
    status='';
    nb=1;
    while nb>0
        nb=COM.BytesAvailable;
        if nb>0
            status=fscanf(COM,'%c',nb);
        end
    end
    fprintf(1,'Activating Shutter: %i of %i\n',n,i)
    fprintf(COM,A);
    status='';
    while numel(strfind(status,'D'))==0
        nb=COM.BytesAvailable;
        if nb>0
            status=fscanf(COM,'%c',nb);
        end
    end
end

disp('shutter closed');
end
