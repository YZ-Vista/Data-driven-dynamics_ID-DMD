clear all, 
clc
load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
currentFolder = pwd;
fullPath = fullfile(currentFolder, 'PIDON');
addpath(fullPath)
load extra_test.mat extra_pred extra_true
Y_test = extra_true.';
Y_pred = extra_pred;
[~,N] = size(Y_test)

for k=1:N
Err(:,k) = abs(Y_pred(:,k)-Y_test(:,k))./max(Y_test(:,k));
end

t = 0:0.01:3;
x = 0:0.01:1; 

subplot(1,3,1)
surf(t(1:3:end),x,Y_test);
shading interp
caxis([-0.5 0.5])
view([0,90])
colormap(subplot(1,3,1),Map_rebu);

subplot(1,3,2)
surf(t(1:3:end),x,Y_pred); 
shading interp
caxis([-0.5 0.5])
view([0,90])
colormap(subplot(1,3,2),Map_rebu);

subplot(1,3,3)
surf(t(1:3:end),x,Err);
shading interp
% caxis([0 6])
view([0,90])
colormap(subplot(1,3,3),Map_grye);
