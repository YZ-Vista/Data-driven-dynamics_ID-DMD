clc
clear
load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
load VORTEX.mat VORTEX PARA CYLIN lx ly;
% VORTALL contains flow fields reshaped into column vectors
X = VORTEX;
[M,N] = size(X{1});
[~,n] = size(X);

r = 120;
EX=[];
A = 3e-1;
for ni=1:N
    randn('state',ni)
    EX(:,ni) = A*randn(1,ly*lx);
end

%% ID-DMD 
start = 40;
for xi=1:n
    Xdn{xi} = X{xi}(:,start:end);
    X{xi} = Xdn{xi}+EX(:,start:end);
end

[Mx,Nx] = size(X{1});
ntrunc = round(1*Nx);

scal = 1;
indx=[1 3 5 7 9 11];
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:ntrunc-1);
    Xc_prime{i} = X{k}(:,2:ntrunc);
    P{i} = PARA{k}*scal;
end

s = 6;
pt = PARA{s}*scal;
run = 1;

for tr = 1:run
    [Phi, Lambda, b, Ubaser] = DMD_for_D1_Bag(Xc,Xc_prime,P,r,pt);
    
    dt = 0.5; % dt can be set as any value as it will not affect the final results
    tspan = [dt:dt:Nx*dt];
    
    XPm(:,1)=Xdn{s}(:,1);
    fhandle = plotCylinderX(real(reshape(XPm(:,1),ly,lx)),1)
    colormap(subplot(1,3,1),Map_rebu)
    axis equal off; drawnow
    b = b\(Ubaser'*Xdn{s}(:,1));

    Err = 0;
    for k = Nx-30 % Time series 2:Nx
        omega = log(Lambda)/dt;
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
        colormap(subplot(1,3,1),Map_rebu)
        axis equal off; drawnow
        
        Vtest = Xdn{s}(:,k);
        fhandle = plotCylinderX(reshape(Vtest,ly,lx),2)
        colormap(subplot(1,3,2),Map_rebu)
        axis equal off; drawnow
        %%%%%%%%%%% err
        Err = Err + abs(Y-Vtest)./max(abs(Vtest));
    end
end
    
AErr = Err/1; % or Err/Nx for all time series
fhandle = plotCylinderX_E(reshape(AErr,ly,lx),3)
colormap(subplot(1,3,3),Map_grye)
axis equal off; drawnow
caxis([0 0.05])
