clear all, 
clc
addpath('D:\Backup\Working papers\2024 DMD for design\GitHub_upload\7_ID_DMD_Droplet\DataSet');

load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
load DropLet_D.mat Xdroplet Vol lx ly;

X = Xdroplet;
PARA = Vol;
r = 200;
[~,n] = size(X);
[Mx,Nx] = size(X{1})
EX=[];
A=0e-5;
for ni=1:Nx
randn('state',ni)
EX(:,ni) = 2*A*randn(1,Mx)-A;
end
for xi=1:n
    Xdn{xi} = X{xi};
    X{xi} = Xdn{xi}+EX;
end

scal = 100;
indx=[2 4 6]; %
% indx=[1:5]; %
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:end-1);
    Xc_prime{i} = X{k}(:,2:end);
    P{i} = PARA{k}/scal;
end
s = 3; % 2
[M,N] = size(X{s});
pt = PARA{s}/scal;
[Phi, Lambda, b, Ubaser] = DMD_for_D1(Xc,Xc_prime,P,r,pt);

%% Animation
anim = VideoWriter('my_animation.avi');
anim.FrameRate = 30;
open(anim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=1/36;
tspan = [0:dt:(N-1)*dt];
omega = log(Lambda)/dt;
wx = abs(diag(omega));

XPm(:,1)=Xdn{s}(:,1);
    fhandle = plotCylinderX(real(reshape(XPm(:,1),lx,ly)),1)
    colormap(subplot(1,2,1),Map_rebu);
    axis equal off; drawnow 
b=b\(Ubaser'*Xdn{s}(:,1));

Err = 0;
for k = 2:N
%% DMD
    Y = zeros(lx*ly,1);
    w = [];
    for i = 1:r
        if real(omega(i,i))>0
            omegai = omega(i,i)-real(omega(i,i));
        else
            omegai = omega(i,i);
        end
        Y = Y+Phi(:,i)*exp(omegai*tspan(k))*b(i);
        w = [w abs(omegai)];
    end  
    fhandle = plotCylinderX(real(reshape(Y,lx,ly)),1)
    colormap(subplot(1,3,1),Map_rebu);
    axis equal off; drawnow 
    
    Vtest = X{s}(:,k-1);
    fhandle = plotCylinderX(reshape(Vtest,lx,ly),2)
    colormap(subplot(1,3,2),Map_rebu);
    axis equal off; drawnow 
    
    frame = getframe(gcf);
    writeVideo(anim, frame);
%%%%%%%%%%%%% err
    Err = Err + abs(Y-Vtest)./max(abs(Vtest));
end
close(anim);

AErr = Err/N;
fhandle = plotCylinderX(reshape(AErr,lx,ly),3)
colormap(subplot(1,3,3),Map_grye);
axis equal off; drawnow
