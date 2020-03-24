test = createEmptyCG2D(0.1, 5, 5);
test.data(25,25) = 1;
result = convFluence(test, 1);

figure;
subplot(1,2,1);
test.plotSlice
colorbar
title('Original');
subplot(1,2,2);
result.plotSlice
colorbar
title('Filtered');