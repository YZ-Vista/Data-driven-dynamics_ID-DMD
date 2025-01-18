clear all, 
clc
addpath('D:\Backup\Working papers\2024 DMD for design\GitHub_upload\6_ID_DMD_Smoke\DataSet');
% 
%% Radius
load smokeD_r1.mat smokeD_r1;
load smokeD_r2.mat smokeD_r2;
load smokeD_r3.mat smokeD_r3;
load smokeD_r4.mat smokeD_r4;
load smokeD_r5.mat smokeD_r5;
load smokeD_r6.mat smokeD_r6;
load smokeD_r7.mat smokeD_r7;
load smokeD_r8.mat smokeD_r8;
load smokeD_r9.mat smokeD_r9;
load smokeD_r10.mat smokeD_r10;
load smokeD_r11.mat smokeD_r11;
% VORTALL contains flow fields reshaped into column vectors
[td,lx,ly] = size(smokeD_r1);

for i = 1:td
    Smoke2D1(:,:) = smokeD_r1(i,:,:);
    Xsmoke{1}(:,i) = reshape(Smoke2D1,lx*ly,1);
    PARA{1} = 5;
    
    Smoke2D2(:,:) = smokeD_r2(i,:,:);
    Xsmoke{2}(:,i) = reshape(Smoke2D2,lx*ly,1);
    PARA{2} = 5.1;
    
    Smoke2D3(:,:) = smokeD_r3(i,:,:);
    Xsmoke{3}(:,i) = reshape(Smoke2D3,lx*ly,1);
    PARA{3} = 5.2;
 
    Smoke2D4(:,:) = smokeD_r4(i,:,:);
    Xsmoke{4}(:,i) = reshape(Smoke2D4,lx*ly,1);
    PARA{4} = 5.3;
    
    Smoke2D5(:,:) = smokeD_r5(i,:,:);
    Xsmoke{5}(:,i) = reshape(Smoke2D5,lx*ly,1);
    PARA{5} = 5.4;
    
    Smoke2D6(:,:) = smokeD_r6(i,:,:);
    Xsmoke{6}(:,i) = reshape(Smoke2D6,lx*ly,1);
    PARA{6} = 5.5;

    Smoke2D7(:,:) = smokeD_r7(i,:,:);
    Xsmoke{7}(:,i) = reshape(Smoke2D7,lx*ly,1);
    PARA{7} = 5.6;
    
    Smoke2D8(:,:) = smokeD_r8(i,:,:);
    Xsmoke{8}(:,i) = reshape(Smoke2D8,lx*ly,1);
    PARA{8} = 5.7;
    
    Smoke2D9(:,:) = smokeD_r9(i,:,:);
    Xsmoke{9}(:,i) = reshape(Smoke2D9,lx*ly,1);
    PARA{9} = 5.8;
    
    Smoke2D10(:,:) = smokeD_r10(i,:,:);
    Xsmoke{10}(:,i) = reshape(Smoke2D10,lx*ly,1);
    PARA{10} = 5.9;
    
    Smoke2D11(:,:) = smokeD_r11(i,:,:);
    Xsmoke{11}(:,i) = reshape(Smoke2D11,lx*ly,1);
    PARA{11} = 6;
end

save Smoke_design_R.mat Xsmoke PARA lx ly

