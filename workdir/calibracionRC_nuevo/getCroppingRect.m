function croppingRect = getCroppingRect(imgPath, msg)
    I=imread(imgPath);
    warning(msg);
    figure(77);
    [~, croppingRect] = imcrop(I);
    close(77);
end

