clear all, 
clc
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

start = 40;
for xi=1:n
    Xdn{xi} = X{xi}(:,start:end);
    X{xi} = Xdn{xi}+EX(:,start:end);
end

for unc = 1:30
indx=[1 3 5 7 9 11];
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:end-1);
    Xc_prime{i} = X{k}(:,2:end);
    P{i} = PARA{k}*1;
end
[Atilde,Btilde,Ubaser] = DMD_for_DPre(Xc,Xc_prime,P,r);

Xdn = X;
tic
Theta_set = -[2:1:12]*pi/180;
for s = 1:length(Theta_set)
    theta = Theta_set(s);
    [M,N] = size(Xdn{s}(:,1:end-1));
    pt = -theta*1;
    CombM = Atilde+pt.*Btilde;
    [W,Lambda] = eig(CombM);         % Step 3
    Phi =  Ubaser*W;       % Step 4
    b = (W*Lambda);

    XPm(:,1) = Xdn{s}(:,1);
    VOE(:,1) = XPm(:,1);
    b = b\(Ubaser'*Xdn{s}(:,1));

    dt=0.5;
    tspan = [dt:dt:N*dt];
%% Generate new data from DMD for design    
  for k=2:N
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

    VOE(:,k) = Y;
  end    
    VORTEX_D{s} = real(VOE);
    PARA_D{s} = -theta;  
end
toc
[Theta, Eness(unc,:)] = Uncertain_Design(VORTEX_D, PARA_D, lx, ly);
EnessRate(unc,:) = (Eness(unc,:));%./Eness(unc,1)
unc
end
% 
% boxplot(EnessRate, Theta);
% hold on
save Ene.mat Eness EnessRate Theta