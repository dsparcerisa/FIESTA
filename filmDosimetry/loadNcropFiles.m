function I = loadNcropFiles(filePaths, nCrops)
% Using imcrop and uiopen dialog box, opens a given number of files
    I = {};
    N = numel(filePaths);
    
    for i=1:N
        %% Load image
        fullImage=imread(filePaths{i});
        
        %% Quick fix for incorrectly saved images
        if size(fullImage,3)>3
            fullImage = fullImage(:,:,1:3);
        end
        
        fprintf('Image %i loaded.\n', i);
        
        %% Crop images
        Npoints=nCrops(i);
        for j=1:Npoints
            figure(1);
            I{(numel(I)+1)}=imcrop(fullImage);
            close(1);
        end
    end
end

