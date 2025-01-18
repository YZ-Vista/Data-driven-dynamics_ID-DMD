clear all, 
clc
load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
currentFolder = pwd;
fullPath = fullfile(currentFolder, 'PIDON');
addpath(fullPath)
load inter_test.mat inter_pred inter_true
Y_test = inter_true.';
Y_pred = inter_pred;
[~,N] = size(Y_test)

for k=1:N
Err(:,k) = abs(Y_pred(:,k)-Y_test(:,k))./max(Y_test(:,k));
end

t = 0:0.01:1;
x = 0:0.01:1; 

subplot(1,3,1)
surf(t,x,Y_test);
shading interp
caxis([-0.5 0.5])
view([0,90])
colormap(subplot(1,3,1),Map_rebu);

subplot(1,3,2)
surf(t,x,Y_pred); 
shading interp
caxis([-0.5 0.5])
view([0,90])
colormap(subplot(1,3,2),Map_rebu);

subplot(1,3,3)
surf(t,x,Err);
shading interp
% caxis([0 0.6])
view([0,90])
colormap(subplot(1,3,3),Map_grye);
