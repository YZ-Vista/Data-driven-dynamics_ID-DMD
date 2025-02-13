clc
clear all
% close all
load BUILDING.mat Y PARA X_T;

n = length(Y{1}(1,:));
for j = 1:6
X{j} = Y{j};
end

[~,nm] = size(X);
for xi = 1:nm
    X{xi} = X{xi}.';
    X_sec{xi} = X{xi}(:,2:3);
    X_T{xi} = X_T{xi}.';
end

indx = [1 3 5 6]; %
term_degree = 8;
ny = 2; % Time delay for output y (k-1,k-2);
maxlag = ny;
for i = 1:length(indx)
    k=indx(i);
    [Xc{i}, Xc_prime{i}] = StateMatricB(X{k},X_sec{k},ny,term_degree);
    P{i} = PARA{k};
end

s = 2; % 2
stat = 3;
[Xt, ~] = StateMatricB(X{s},X_sec{s},ny,term_degree);
pt = PARA{s};
%% DMD for design
% [m,n] = size(Xt);
r = 100
[Phi, Lambda, b, Ubaser, Mcom] = DMD_for_D1(Xc,Xc_prime,P,r,pt);
b=b\(Ubaser'*Xt(:,1));

%% Plot time series
fs = 64; dt = 1/fs;
figure
plot((1:n-4)*dt,Y{s}(stat,3:end-maxlag))
hold on

%% Mode prediction
omega = log(Lambda)/dt;
tspan = [0:1/fs:(n-1)/fs];
Ydmd = zeros(1,n);
w = [];
for i=1:r
    if real(omega(i,i))>0
        omegai = omega(i,i)-real(omega(i,i));
    else
        omegai = omega(i,i);
    end
    Ydmd = Ydmd+Phi(2*stat-1,i)*exp(omegai*tspan)*b(i);
    w = [w imag(omegai)];
end
hold on
plot((1:n-4)*dt,Ydmd(2:end-maxlag-1))

Err = abs(Y{s}(stat,3:end-maxlag)-Ydmd(2:end-maxlag-1))./max(abs(Y{s}(stat,3:end-maxlag,1)));
AErr = sum(Err)/n
