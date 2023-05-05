
function xlevel = findEdge(x,y,level)

    [~, maxPos] = max(y);
    yNorm = y / max(y);
    yNorm = yNorm(1:maxPos);
    xNorm = x(1:maxPos);
    
    % find first usable point
    firstUsablePoint = 1 + find(diff(yNorm(1:maxPos))<0,1,'last');
    if isempty(firstUsablePoint)
        firstUsablePoint = 1;
    end
    intVal = firstUsablePoint:maxPos;
    xlevel = interp1(yNorm(intVal),x(intVal),level);
   
end