function [status, posX, posY, posZ] = monitorStatus(COM, pauseTime)
% 0: off, 1: ready, 2: moving

if ~exist('pauseTime')
    pauseTime = 0.1;
end
step2cm = 2.5 / 10000;
posX = nan;
posY = nan;
posZ = nan;

statusLetter = readStatus(COM,'V',pauseTime);

switch statusLetter
    case 'R'
        status = 1;
        while isnan(posX)
            posX = step2cm * str2double(readStatus(COM,'X',pauseTime));
        end
        while isnan(posY)        
            posY = step2cm * str2double(readStatus(COM,'Y',pauseTime));
        end
        while isnan(posZ)        
            posZ = step2cm * str2double(readStatus(COM,'Z',pauseTime));  
        end
    case 'B'
        status = 2;
        while isnan(posX)
            posX = step2cm * str2double(readStatus(COM,'X',pauseTime));
        end
        while isnan(posY)
            posY = step2cm * str2double(readStatus(COM,'Y',pauseTime));       
        end
    otherwise
        [status, posX, posY, posZ] = monitorStatus(COM, pauseTime);
end

end

