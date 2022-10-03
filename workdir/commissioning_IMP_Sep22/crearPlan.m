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
Xpositions = -5:0.2:5;
Ypositions = -3.6:0.2:3.6;

%% We now want to interpolate our image to the resolution that we can achieve with a spacing of 0.1 cm, 71x81 points.
logoResize = imresize(((logoImBw)),[numel(Ypositions) numel(Xpositions)]);
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
noZeroMask = logoResize(:,:,1)>15;
figure; imagesc(noZeroMask); title('Non-zero mask'); colormap('gray');

%% And finally, we copy those values into the Plan structure
planName = 'plan0.txt';
planPath = fullfile(Fcfg.planPath, planName);
thePlan = readPlan(planPath);

thePlan.numSpots = numel(newBeamPosX(noZeroMask));
thePlan.Z = newBeamPosX(noZeroMask)';
thePlan.Y = -newBeamPosY(noZeroMask)';
thePlan.X = 0 * thePlan.Z;
thePlan.T_out_us = 300 * ones(size(thePlan.X));
thePlan.Nshots = round(logoResize(noZeroMask)/10);

planNameGFN = 'planGFN.txt';
planPathGFN = fullfile(Fcfg.planPath, planNameGFN);
writePlan(thePlan,planPathGFN)
%% 