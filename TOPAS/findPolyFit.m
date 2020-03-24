function pval = findPolyFit(Zvalues,Svalues)

%% Create fit and plot
F2 = fit((ZValues(minValidIndex:maxValidIndex))', (Sflu(minValidIndex:maxValidIndex))', 'poly2', 'Lower', [0 0 0], 'Upper', [Inf Inf 0]);
pflu = [F2.p1 F2.p2 F2.p3];

figure
plot(ZValues(minValidIndex:maxValidIndex), Sflu(minValidIndex:maxValidIndex), 'b.');
hold on
plot(ZValues(1:maxValidIndex), polyval(pflu, ZValues(1:maxValidIndex)), 'r-');
axis([ZValues(1) ZValues(maxValidIndex) 0 polyval(pflu, ZValues(maxValidIndex))*1.1]);
%%
% maxFitIgnored = maxFitIgnored +NZfluMin;
% Sflu(1:maxFitIgnored) = dR;
% fprintf('Polynomial Fluence: Fittingin the range %4.2f <= Z < %4.2f\n', Zvaluesflu(maxFitIgnored), Zvaluesflu(length(Zvaluesflu)));
% 
% % pflu=polyfit(Zvalues,Sflu,2); %% polyfit no permite restringir valores
% F = fit(Zvaluesflu(maxFitIgnored+1:length(Zvaluesflu))', Sflu(maxFitIgnored+1:length(Sflu))', 'poly2', 'Lower', [0 0 0]);
% pflu = coeffvalues(F);



end

