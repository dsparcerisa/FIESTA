function [steps] = linearstage(COM,id,direction,distance)
    mm2steps=1000/2.5;
    if direction==1 
        strdir='';
    elseif direction==-1
        strdir='-';
    elseif direction==0
        strdir='';
        if distance~=0
            distance
            display('wrong distance for 0 direction!')
            pause
        end
    else
        direction
        display('wrong direction!')
        pause
    end
    steps=round(distance*mm2steps);
    if steps>0
%    fprintf(1,'Moving stage #%d: %d steps %f mm, direction: %d\n',id,steps,distance,direction);
    command=strcat('C,I',num2str(id),'M',strdir,num2str(steps),',R');%Set steps to move (2.5 mm / 1000 steps)
    fprintf(COM,command);
    %fprintf(COM,'C,V,R'); 
%    status='';
%     disp('Reading status of linear stage');
%     while numel(strfind(status,'^'))==0
%         nb=COM.BytesAvailable;
%         if nb>0
%             status=fscanf(COM,'%c',nb);
%         end
%     end

%    status = readStatus(COM);
%    while numel(status)==0
%        status = readStatus(COM)
%        disp(status)
%    end
    
%    disp('done moving linear stage');
%     if id==1
%         savezpos(direction,steps);
%     elseif id==2 
%         savexpos(direction,steps);
%     elseif id==3 
%         saveypos(direction,steps);
%     else
%         id
%         disp('Wrong stage ID');
%     end
    end
end

