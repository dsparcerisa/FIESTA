function limits = autoCenter(COM)
    limits = nan(1,6);
    
    oldLimits = findLimits(COM);
    
    centerX = 0.5*(oldLimits(2)+oldLimits(1));
    limits(1) = oldLimits(1) - centerX;
    limits(2) = oldLimits(2) - centerX;
    
    centerY = 0.5*(oldLimits(4)+oldLimits(3));
    limits(3) = oldLimits(3) - centerY;
    limits(4) = oldLimits(4) - centerY;   

    centerZ = 0.5*(oldLimits(6)+oldLimits(5));
    limits(5) = oldLimits(5) - centerZ;
    limits(6) = oldLimits(6) - centerZ;

    % Round to the nearest mm
    limits = round(limits,1);
    
    while(true)
        status = monitorStatus(COM);
        if status==1
            break
        end
    end
    
    finished = stageControl_moveToAbsPos(COM, [centerX centerY centerZ])
    while(finished == false)
        if finished
            break
        end
    end    
    
    % Zero here
    readStatus(COM, 'N');
    
    % Final stat~us verification
    [~, posX, posY, posZ] = monitorStatus(COM);
    fprintf('Zeroed. Current position is (%3.2f, %3.2f, %3.2f)\n', posX, posY, posZ);
end

