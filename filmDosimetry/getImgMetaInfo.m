function [pixelsXcm, maxInt] = getImgMetaInfo(filePath)
    imgInfo = imfinfo(filePath);
    % 1 dpi = 0.393701 pixel/cm
    pixelsXcm = imgInfo.XResolution * 0.393701;
    maxInt = imgInfo.MaxSampleValue(1);
end

