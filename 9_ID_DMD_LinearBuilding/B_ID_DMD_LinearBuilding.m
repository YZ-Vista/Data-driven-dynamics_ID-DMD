clc
clear
load BUILDING.mat BUILDING PARA X_T;
X = BUILDING;
r = 16 
indx=[1 3 5 6];
[m,n] = size(X{1});
ntrunc = round(1*n);
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:ntrunc-1);
    Xc_prime{i} = X{k}(:,2:ntrunc);
    P{i} = PARA{k};
end
s = 4;
pt = PARA{s};
[Phi, Lambda, b] = DMD_for_D1(Xc,Xc_prime,P,r,pt);
% 
%% Code to plot the modes 
for i=1:8
    figure
    PlotMode(real(Phi(1:4:end,i)));
    hold on
    PlotMode(real(Phi(2:4:end,i)));
    hold on
    PlotMode(real(Phi(3:4:end,i)));
    hold on
    PlotMode(real(Phi(4:4:end,i)));
end

fs=80; dt=1/fs; np=100; N=np*fs;
figure
plot((1:n-1)*dt,X_T{s}(1,5:end))
hold on

%% Mode prediction
omega = log(Lambda)/dt;
tspan = [0:1/fs:(N-1)/fs];
Y = zeros(1,N);
w = [];
for i=1:r
    if real(omega(i,i))>0
        omegai = omega(i,i)-real(omega(i,i));
    else
        omegai = omega(i,i);
    end
    Y = Y+Phi(1,i)*exp(omegai*tspan)*b(i);
    w = [w abs(omegai)];
end
hold on
plot((1:n-1)*dt,Y(2:end-2))

Err = abs(X_T{s}(1,5:end)-Y(2:end-2))./max(abs(X_T{s}(1,5:end)));
AErr = sum(Err)/n

% % hold on
