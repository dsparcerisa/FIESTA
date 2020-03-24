function beamProfile = createPinHole(dxy, sizeX, sizeY, pinholeRadius)
% CG2D beamProfile = createPinHole(double dxy, double sizeX, 
%   double sizeY, double pinholeRadius) 
% Creates a CG2D object contaning a beam profile with open/closed pixels
% simulating a pinhole collimator. Sum is NORMALIZED to 1.
beamProfile = createEmptyCG2D(dxy, sizeX, sizeY);
Xval = beamProfile.getAxisValues('X');
Yval = beamProfile.getAxisValues('Y');
distanceMatrix = sqrt(Xval.^2 + (Yval').^2);
beamProfile.data = distanceMatrix<=pinholeRadius;
beamProfile.data = beamProfile.data / sum(beamProfile.data(:));
end

