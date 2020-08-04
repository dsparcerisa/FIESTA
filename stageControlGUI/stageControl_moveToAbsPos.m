function [finished, finalPos] = stageControl_moveToAbsPos(COM, absPos)
% DOES NOT VERIFY absPos is feasible

finished = false;
pauseTime = 0.05;
% Move stage to position
[status, posX, posY, posZ] = monitorStatus(COM, pauseTime);

if status ~=1
    error('Stage is not in READY status');
end

deltaX_mm = 10*(absPos(1) - posX);
deltaY_mm = 10*(absPos(2) - posY);
deltaZ_mm = 10*(absPos(3) - posZ);

if deltaX_mm ~= 0
    %tic
    linearstage(COM,1,sign(deltaX_mm),abs(deltaX_mm));  
    while(true)
        [status, posX, posY, posZ] = monitorStatus(COM, pauseTime);
        if status==1
            break
        end
    end
    %toc
end

if deltaY_mm ~= 0
    %tic;
    linearstage(COM,2,sign(deltaY_mm),abs(deltaY_mm));    
    while(true)
        [status, posX, posY, posZ] = monitorStatus(COM, pauseTime);
        if status==1
            break
        end
    end
    %toc;
end

if deltaZ_mm ~= 0
    %tic;
    linearstage(COM,3,sign(deltaZ_mm),abs(deltaZ_mm));    
    while(true)
        [status, posX, posY, posZ] = monitorStatus(COM, pauseTime);
        if status==1
            break
        end
    end
    %toc;
end

finished = true;
finalPos = [posX, posY, posZ];
if (absPos(1) ~= posX)
    warning('Could not reach desired X position: wanted %3.3f, arrived at %3.3f', absPos(1), posX)
end
if (absPos(2) ~= posY)
    warning('Could not reach desired Y position: wanted %3.3f, arrived at %3.3f', absPos(2), posY)
end
if (absPos(3) ~= posZ)
    warning('Could not reach desired Z position: wanted %3.3f, arrived at %3.3f', absPos(3), posZ)
end
end

