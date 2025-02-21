load Bur_Test_Extra.mat x t Burg_pre Burg_true pt
parameter = pt

subplot(1,2,1)
surf(t(1:3:end),x,Burg_pre(:,1:3:end)), 
% shading interp, lighting phong, axis tight
% colormap(summer); set(gca,'zlim')
% light('color',[1 1 0],'position',[-1,2,2])
% material([0.30 0.60 0.60 40.00 1.00]);
% title('The Burgers equation')
% xlabel('t'); ylabel('x')
% view([0,90])
view([90,-90])
colorbar

subplot(1,2,2)
surf(t(1:3:end),x,Burg_true(:,1:3:end)), 
% shading interp, lighting phong, axis tight
% colormap(summer); set(gca,'zlim')
% light('color',[1 1 0],'position',[-1,2,2])
% material([0.30 0.60 0.60 40.00 1.00]);
% title('The Burgers equation')
% xlabel('t'); ylabel('x')
% view([0,90])
view([90,-90])
colorbar

Err= abs(Burg_pre-Burg_true);

figure
surf(t,x,Err), 
% shading interp, lighting phong, axis tight
% colormap(summer); set(gca,'zlim')
% light('color',[1 1 0],'position',[-1,2,2])
% material([0.30 0.60 0.60 40.00 1.00]);
% title('The Burgers equation')
% xlabel('t'); ylabel('x')
view([90,-90])
colorbar
