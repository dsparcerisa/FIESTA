function motorResponse = readStatus(s2, in, pauseTime)
    if ~exist('pauseTime')
        pauseTime = 0.05;
    end
    if exist('in')
        % disp(in)
        fprintf(s2, in);
        % pause(0.05);        optimal
        pause(pauseTime)
    end
    if s2.BytesAvailable>0
        motorResponse = fscanf(s2,'%c',s2.BytesAvailable);
        motorResponse = motorResponse(motorResponse~='^');
    else
        motorResponse = '';
    end
end

