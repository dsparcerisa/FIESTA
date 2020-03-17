function startAutoMonitorStatus(app, pauseTime)

port = app.getPort;

while(true)
    [status, x, y, z] = monitorStatus(port);
    app.setStatus(status);
    switch status
        case 1
            app.updatePosition([x y z]);

        case 2
            oldPos = app.getCurrentPos;
            app.updatePosition([ x y oldPos(3)]);
    end
    pause(pauseTime);
end

end

