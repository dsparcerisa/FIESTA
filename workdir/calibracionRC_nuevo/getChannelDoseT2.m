function D = getChannelDoseT2(CoefX2, PVX)
    % CoefX2 = [p1 p2 q1]
    p1 = CoefX2(1);
    p2 = CoefX2(2);
    q1 = CoefX2(3);
    ODX = -log10(PVX);
    D = (p2 - ODX.*q1) ./ (ODX - p1);
end

