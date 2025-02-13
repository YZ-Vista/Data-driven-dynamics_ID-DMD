clc
clear
load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
load Cavity_D.mat Y_Cavity PARA Nx Ny;
lx = Nx; ly = Ny
X = Y_Cavity;
r = 60;
[~,nm] = size(X);
[M,N] = size(X{1});
ntrunc = round(1*N);
EX=[];
A=0e-2;
for ni=1:N
randn('state',ni)
EX(:,ni) = 2*A*randn(1,M)-A;
end
for xi=1:nm
    Xdn{xi} = X{xi};
    X{xi} = Xdn{xi}+EX;
end

indx = [1:4:25]; % indx = [1:2:16];
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:ntrunc-1);
    Xc_prime{i} = X{k}(:,2:ntrunc);
    P{i} = PARA{k};
end
s = 7; % 24
pt = PARA{s};
[Phi, Lambda, b, Ubaser, Mcom] = DMD_for_D1(Xc,Xc_prime,P,r,pt);
%% Animation
anim = VideoWriter('CavityFlow_animation.avi');
anim.FrameRate = 30;
open(anim);

dt=0.01;
tspan = [0:dt:(N-1)*dt];
omega = log(Lambda)/dt;
wx = abs(diag(omega));

XPm(:,1)=Xdn{s}(:,1);
    fhandle = plotCylinderX(real(reshape(XPm(:,1),ly,lx)),2,1)
    colormap(subplot(1,2,1),Map_rebu);
    axis equal off; drawnow 
b=b\(Ubaser'*X{s}(:,1));

%% DMD
Err = 0;
for k=2:N %
    Y = zeros(ly*lx,1);
    w = [];
    for i=1:r
        if real(omega(i,i))>0
            omegai = omega(i,i)-real(omega(i,i));
        else
            omegai = omega(i,i);
        end
        Y = Y+Phi(:,i)*exp(omegai*tspan(k))*b(i);
        w = [w abs(omegai)];
    end  
    fhandle = plotCylinderX(real(reshape(Y,ly,lx)),2,1)
    colormap(subplot(1,2,1),Map_rebu);
    axis equal off; drawnow 
    
    Vtest = X{s}(:,k-1);
    fhandle = plotCylinderX(reshape(Vtest,ly,lx),2,2)
    colormap(subplot(1,2,2),Map_rebu);
    axis equal off; drawnow 

    frame = getframe(gcf);
    writeVideo(anim, frame);
end
close(anim);
