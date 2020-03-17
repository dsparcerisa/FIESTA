function stageControl_setPosition(app, position)
% moves stage to defined position
if app.isWithinLimits(position)
    message = sprintf('Move to position -> (X: %3.2f  Y: %3.2f  Z: %3.2f)\n', position(1), position(2), position(3));
    app.logLine(message)
           
    % Move stage to position
    currentPos = app.getCurrentPos;
    deltaX_mm = 10*(position(1) - currentPos(1));
    deltaY_mm = 10*(position(2) - currentPos(2));
    deltaZ_mm = 10*(position(3) - currentPos(3));
    if deltaX_mm ~= 0
        %tic
        linearstage(app.getPort,1,sign(deltaX_mm),abs(deltaX_mm));
        %toc
    end
    if deltaY_mm ~= 0        
        %tic;
        linearstage(app.getPort,2,sign(deltaY_mm),abs(deltaY_mm));
        %toc;
    end
    if deltaZ_mm ~= 0        
        %tic;
        linearstage(app.getPort,3,sign(deltaZ_mm),abs(deltaZ_mm));
        %toc;
    end
else
    message = 'ERROR: Target position is not within defined limits';
    app.logLine(message);
end
end


