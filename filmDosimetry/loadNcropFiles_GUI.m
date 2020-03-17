function I = loadNcropFiles_GUI(N, Fcfg)
% Using imcrop and uiopen dialog box, opens a given number of files
    I = {};
    for i=1:N
        %% Load image
        pathModel = [Fcfg.filmPath filesep '*.tif'];
        [fileName, filePath] = uigetfile(pathModel);
        imgPath = [filePath, fileName];
        fullImage=imread(imgPath);
        
        %% Quick fix for incorrectly saved images
        if size(fullImage,3)>3
            fullImage = fullImage(:,:,1:3);
        end
        
        fprintf('Image %i loaded.\n', i);
        
        %% Crop images
        Npoints=input('Number of cropped images for this image: ');
        for j=1:Npoints
            figure(1);
            I{(numel(I)+1)}=imcrop(fullImage);
            close(1);
        end
    end
end

