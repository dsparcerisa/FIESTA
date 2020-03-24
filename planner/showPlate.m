function showPlate(plate)
figure
NX = size(plate, 1);
NY = size(plate, 2);
[x,y] = meshgrid(1:NY, 1:NX);
scatter(x(:), y(:), 300, plate(:), 'filled')
set(gca, 'Ydir', 'reverse')
hold on
scatter(x(:), y(:), 600)
axis([0 (NY+1) 0 (NX+1)]);
colorbar
end

