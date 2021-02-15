function croppedI = cropCentralStrip(I,wdPx,htPx)
%Crop central strip of wdPx x htPx pixels of image I

totalHt = size(I,1);
totalWd = size(I,2);
X_minus = totalWd/2 - wdPx/2;
Y_minus = totalHt/2 - htPx/2;

croppedI = imcrop(I, [X_minus Y_minus wdPx htPx]);

end

