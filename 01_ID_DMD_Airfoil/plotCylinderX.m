function fhandle = plotCylinderX(VORT,i)

fhandle = figure(1)
subplot(1,2,i)

contourf(VORT,-20:0.05:20, 'LineColor','none')
% clim([-19 4])
caxis([-2 2]) % clim after MATLAB 2022
colorbar;
