function D = getChannelDoseT1(CoefX, PVX)
    D = CoefX(3) + CoefX(2) ./ (PVX - CoefX(1));
end

