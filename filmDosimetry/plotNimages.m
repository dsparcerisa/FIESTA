function plotNimages(I)
    H = floor(sqrt(numel(I)));
    W = ceil(numel(I) / H);
    for i=1:numel(I)
        subplot(H, W, i);
        if(numel(size(I{i})) == 3)
            imshow(I{i});
        else
            imagesc(I{i});
                caxis([0 10]);
        end
        title(sprintf('%i',i));
    end
end

