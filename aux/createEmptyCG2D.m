function emptyCG2D = createEmptyCG2D(dxy, sizeX, sizeY)
% CG2D emptyCG2D = createEmptyCG2D(double dxy, double sizeX, double sizeY)
NX = round(sizeX/dxy);
NY = round(sizeY/dxy);
maxX = (NX*dxy/2) - dxy/2;
maxY = (NY*dxy/2) - dxy/2;
emptyCG2D = CartesianGrid2D(-maxX,maxX,-maxY,maxY,dxy,dxy,NX,NY);
end

