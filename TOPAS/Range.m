function NZ=Range(E)
%Calculation of the range for an array of Edep along the Z axis
NZ=0;
for i=1:length(E);    
    if E(i)==0;
        break
    end
    NZ=NZ+1;
end
    