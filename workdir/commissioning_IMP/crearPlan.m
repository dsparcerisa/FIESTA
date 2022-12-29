%% Now we will load an image and convert it to greyscale:
%logoIm = imread(fullfile(testPWD,'aux_files','LOGO.png'));
logoIm = imread('GFNsimplif.png');
%logoIm = imread('silueta.jpg');
logoImBw = logoIm; % 255 - min(logoIm,[],3);
imagesc(logoImBw)
colormap('gray');
title('Example image');

%% Our spots in the calculation grid will be placed in these positions
% X direction (BEV) --> From -4 to 4 in intervals of 0.1 cm
% Y direction (BEV) --> From -3.5 to 3.5 in intervals of 0.1 cm
Xpositions = -4:0.1:4;
Ypositions = -3.5:0.1:3.5;

%% We now want to interpolate our image to the resolution that we can achieve with a spacing of 0.1 cm, 71x81 points.
logoResize = imresize(flipud(rot90(logoImBw)),[numel(Ypositions) numel(Xpositions)]);
imagesc(Xpositions,Ypositions,logoResize)
colormap('gray')
title('Rotated and interpolated image');
xlabel('X positions');
ylabel('Y positions');

% (the rotation is performed so that, after applying the nominal couch and
% gantry rotations, the final dose distribution has the right orientation.)

%% The key thing is to modify the current spot list and produce a new spot list with the following characteristics:
% 1. X and Y positions have to match our calculation grid (7x8 cm with 0.1 cm spacing)
% 2. The beam weights must be proportional to the image intensity
% 3. Pixels with zero intensity must be removed from the spot list.

%% We use the function meshgrid (type help meshgrid if you are unfamiliar with it) to create coordinate points for X and Y positions
[newBeamPosX, newBeamPosY] = meshgrid(Xpositions, Ypositions);

%% We create a mask to identify non-zero points, so that we use only those spots
noZeroMask = logoResize>0;
figure; imagesc(noZeroMask); title('Non-zero mask'); colormap('gray');

%% And finally, we copy those values into the Plan structure
FM.plan.fields(1).maxSpots = numel(newBeamPosX(noZeroMask));
FM.plan.fields(1).beamPosX = newBeamPosX(noZeroMask)';
FM.plan.fields(1).beamPosY = newBeamPosY(noZeroMask)';
FM.plan.fields(1).beamWeights = double(50+logoResize(noZeroMask))';
FM.plan.fields(1).generateBeamIndexes(FM.DG);
figure; FM.plan.fields(1).plotSliceSpots(1); colorbar;
FM.plan.showPlan

%% Before starting the calculation, we select a spot size of 1mm FWHM
FM.plan.fields(1).spotSizesX = 0.1;
FM.plan.fields(1).spotSizesY = 0.1;

% This number is actually unrealistic (clinical spot sizes are the order of
% 0.5-1cm in air FWHM), but it will produce a much sharper image in our calculation.



%% This is what the initial dose distribution looks like
figure; subplot(1,2,1); 
FM.plotDoseSlice('Y',2.5)
axis image
title('Transversal dose distribution');
subplot(1,2,2)
axis image
FM.plotDoseSlice('X',0);
title('Longitudinal dose distribution');

%% Tasks and questions:

% 1. Find an image with low resolution and approximately square.
% Replace the FoCa logo and re-run the dose calculation with your own
% image.

% 2. Create plots at depths = 1cm (X=-6cm), 9.25 cm (X=2.25) and 9.5 cm (X=2.5 cm). Which one has a better
% resolution? How do these different resolutions relate to the position of
% the peak? What is the projected spot size at each of these depths?
% Hint: Use this code to see the depth dose distribution

% FM.PPD.plotDDD;
% axis([0 12 0 3]);

% 3. Modify the spot sizes to 0.6 cm FWHM in both directions and re-run the
% dose calculation.
% 3a) At which position would you put a radiochromic film in
% order to register the image as clearly as possible?
% 3b) If the total dose is scaled so that it delivers a max dose of 2 Gy
% near the peak position, what will be the maximum dose received by the
% film?
% (Hint: total dose is saved in FM.results.fullDoseGrid )
% (Hint: you can use function FM.results.getSlice(...), type 'help ResultGrids' for help )
% 3c) If the sensitivity of the film is 0.2 Gy (i.e. doses under this
% threshold will not be registered), what will be the resulting image?