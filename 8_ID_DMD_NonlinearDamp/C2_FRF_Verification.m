clc
clear
np = 500; fs = 32; dt = 1/fs;
N = fs*np;
tspan=[0:1/fs:(N-1)/fs];

k1 = 1e2; c3_set = [1:1:20];
omega = 10;
x0=[0 0];

%% Spectrum
c3 = c3_set(15);
[t,x]=ode45(@odetestT,tspan,x0,[],omega,k1,c3);
y=x(:,1);

Nsamp=fs./(omega/(2*pi))*60;
stat = 1500;
[ftranY,ome,~]=FFTS(y(stat:round(stat+Nsamp)),omega,fs,1,'P');

figure(1)
plot(ome,abs(ftranY))
hold on

figure(2)
plot(t(stat:round(stat+Nsamp)),y(stat:round(stat+Nsamp)))
