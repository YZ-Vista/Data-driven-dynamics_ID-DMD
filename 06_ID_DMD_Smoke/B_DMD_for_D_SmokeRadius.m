clear all, 
clc
load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
load Smoke_design_R.mat Xsmoke PARA lx ly;
% VORTALL contains flow fields reshaped into column vectors
X = Xsmoke;
r = 350; %350;
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
start = 20;
scal = 10;
indx=[1 2 5 6]; %
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,start:end-1);
    Xc_prime{i} = X{k}(:,start+1:end);
    P{i} = PARA{k}/scal;
end
s = 3; % 2
[M,N] = size(X{s});
pt = PARA{s}/scal;
[Phi, Lambda, b, Ubaser, Mcom] = DMD_for_D1(Xc,Xc_prime,P,r,pt);

dt = 1; % arbitrary positive value
tspan = [0:dt:(N-start-1)*dt];
omega = log(Lambda)/dt;
wx = abs(diag(omega));

XPm(:,1)=Xdn{s}(:,start);
    fhandle = plotCylinderX_Denoise(real(reshape(XPm(:,1),ly,lx)),3,1)
    colormap(subplot(1,3,1),Map_rebu);
    axis equal off; drawnow 
b=b\(Ubaser'*Xdn{s}(:,start));
    
%% DMD
Err = 0;
for k = 2:N-start %
    Y = zeros(ly*lx,1);
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
    fhandle = plotCylinderX_Denoise(real(reshape(Y,ly,lx)),3,1)
    colormap(subplot(1,3,1),Map_rebu);
    axis equal off; drawnow 
    
    Vtest = X{s}(:,start+k-1);
    fhandle = plotCylinderX(reshape(Vtest,ly,lx),3,2)
    colormap(subplot(1,3,2),Map_rebu);
    axis equal off; drawnow 

%%%%%%%%%%%%% err
    Err = Err + abs(Y-Vtest)./max(abs(Vtest));
end

AErr = Err/(N-start);
fhandle = plotCylinderX(reshape(AErr,lx,ly),3,3)
colormap(subplot(1,3,3),Map_grye);
axis equal off; drawnow
