function wellMap = getWell(canvas, wellR, wellPos)

wellMap = CartesianGrid2D(canvas);
Xval = (wellMap.getAxisValues('X'))';
Yval = (wellMap.getAxisValues('Y'));
distance2 = (Xval-wellPos(1)).^2 + (Yval-wellPos(2)).^2;
wellMap.data = distance2 <= (wellR^2);

end

