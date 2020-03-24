function showCubeta(plate)
figure
[x,y] = meshgrid(1:2, 1:4);
scatter(x(:), y(:), 300, plate(:), 'filled')
set(gca, 'Ydir', 'reverse')
hold on
scatter(x(:), y(:), 600)
axis([0 3 0 5]);
colorbar
end

