clc
clear
fs=80; dt=1/fs; np=100; N=np*fs;
%% parameters % 1-4: bottom to top
mC1=5*10^6; mC2=4*10^6; mA2=3*10^6; mA1=2*10^6;
kD=[1:0.5:5]*10^9; kC2=2000*10^6; kB=3000*10^6; kA1=1000*10^6;
cC1=1e6; cC2=1e5; cB=1e5; cA1=1e5;
for i=1:6
kC1=kD(i); 
%% Numerical-Random 
tspan=[0:1/fs:N/fs];
x0=zeros(1,8);
x0(1)=1;
[t,x]=ode45(@ode4dof,tspan,x0,[],tspan,...
    cC1,cC2,cB,cA1,kC1,kC2,kB,kA1);

BUILDING{i}=[x(4:end,1) x(3:end-1,1) x(2:end-2,1) x(1:end-3,1)...
    x(4:end,2) x(3:end-1,2) x(2:end-2,2) x(1:end-3,2)...
    x(4:end,3) x(3:end-1,3) x(2:end-2,3) x(1:end-3,3)...
    x(4:end,4) x(3:end-1,4) x(2:end-2,4) x(1:end-3,4)].';
PARA{i} = kC1/1e9;
X_T{i}=x.';
end
save BUILDING.mat BUILDING PARA X_T



