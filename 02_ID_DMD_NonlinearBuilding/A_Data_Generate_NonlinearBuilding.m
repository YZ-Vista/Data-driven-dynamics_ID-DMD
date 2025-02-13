clc
clear

fs=64; dt=1/fs; np=150; N=np*fs;

%% parameters % 1-4: bottom to top
mC1=5*10^6; mC2=4*10^6; mA2=3*10^6; mA1=2*10^6;
kC1=1500*10^6; kC2=2000*10^6; kB=3000*10^6; kA1=1000*10^6;
cC1=1e5; cC2=1e5; cB=1e5; cA1=1e5;
cn3D=[0.1 0.5 0.8 1 3 5]*10^3;

for i=1:6
cn3 = cn3D(i); 
%% Numerical-Random 
tspan=[0:1/fs:N/fs];
x0=zeros(1,8);
x0(1)=1;
[t,x]=ode45(@NB_ode4dof,tspan,x0,[],tspan,...
    cn3,cC1,cC2,cB,cA1,kC1,kC2,kB,kA1);
figure
plot(t,x(:,1))
hold on

plot(t,x(:,2)+2)
hold on

plot(t,x(:,3)+4)
hold on

plot(t,x(:,4)+6)

Y{i}=[x(:,1) x(:,2) x(:,3) x(:,4)].';
PARA{i} = cn3/1e3;
X_T{i}=x.';
end
save BUILDING.mat Y PARA X_T


