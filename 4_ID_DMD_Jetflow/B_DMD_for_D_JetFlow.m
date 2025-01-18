clc
clear
load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
load Jet_D.mat Y_Jet PARA Nx Ny;
lx = Nx; ly = Ny
X = Y_Jet;
r = 150; %150;
[~,nm] = size(X);
[M,N] = size(X{1});
ntrunc = round(N);
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

indx = [1 3 5 7]; % indx = [4 5 6 7];
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:ntrunc-1);
    Xc_prime{i} = X{k}(:,2:ntrunc);
    P{i} = PARA{k};
end
s = 2; % 24
pt = PARA{s};
[Phi, Lambda, b, Ubaser] = DMD_for_D1(Xc,Xc_prime,P,r,pt);

dt=1/36;
tspan = [0:dt:(N-1)*dt];
omega = log(Lambda)/dt;
wx = abs(diag(omega));

XPm(:,1)=Xdn{s}(:,1);
    fhandle = plotCylinderX(real(reshape(XPm(:,1),ly,lx)),1)
    colormap(subplot(1,2,1),Map_rebu);
    axis equal off; drawnow 
b=b\(Ubaser'*X{s}(:,1));

Err = 0;
for k=2:N %-0 -40 -60
    
%% DMD
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
    fhandle = plotCylinderX(real(reshape(Y,ly,lx)),1)
    colormap(subplot(1,3,1),Map_rebu);
    axis equal off; drawnow 
    
    Vtest = X{s}(:,k-1);
    fhandle = plotCylinderX(reshape(Vtest,ly,lx),2)
    colormap(subplot(1,3,2),Map_rebu);
    axis equal off; drawnow 
%%%%%%%%%%%%%% err
    Err = Err + abs(Y-Vtest)./max(abs(Vtest));
end

AErr = Err/N;
fhandle = plotCylinderX(reshape(AErr,lx,ly),3)
colormap(subplot(1,3,3),Map_grye);
axis equal off; drawnow