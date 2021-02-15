function outD = getChannelDoseT3(CoefX3, PVX, PV0X)

    a = CoefX3(1);
    b = CoefX3(2);
    c = CoefX3(3);    
    
    for i=1:numel(PVX)
        netOD(i) = log10(PV0X./PVX(i));
        syms D
        try
            outD(i) = real(double(vpasolve(a*D + b*D^c == netOD(i), D, 1)));
        catch
            outD(i) = nan;
        end
    end
end

