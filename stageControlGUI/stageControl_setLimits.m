function status = stageControl_setLimits(app, limitVector)
% Sets current limits [Xmin Xmax Ymin Ymax Zmin Zmax] in stageControl app
    app.updateLimits(limitVector);
    message = sprintf('Limits -> X: [%3.2f to %3.2f]  Y: [%3.2f to %3.2f]  Z: [%3.2f to %3.2f]\n', limitVector(1), ...
        limitVector(2), limitVector(3), limitVector(4), limitVector(5), limitVector(6));
    app.logLine(message);
    status = true;
end

