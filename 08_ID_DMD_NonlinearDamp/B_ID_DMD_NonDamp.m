clc
clear
load Duff_D.mat Y PARA;

n = length(Y{1});
ntrunc = round(1*n);

for j = 1:5
    X{j} = Y{j}(1:ntrunc);
end

indx = [1 3 5]; %
term_degree = 8;
ny = 2; % Time delay for output y (k-1,k-2);
maxlag = ny;
for i = 1:length(indx)
    k=indx(i);
    [Xc{i}, Xc_prime{i}] = StateMatric(X{k},ny,term_degree);
    P{i} = PARA{k}/10;
end

s = 4; % 2
[Xt, ~] = StateMatric(X{s},ny,term_degree);
pt = PARA{s}/10;
%% DMD for design
r = 30
[Phi, Lambda, b, Ubaser] = DMD_for_D1(Xc,Xc_prime,P,r,pt);

%% Plot time series
fs = 32; dt = 1/fs;
figure
plot((1:n-4)*dt,Y{s}(3:end-maxlag))
hold on

XPm(:,1)=Xt(:,1);
b=b\(Ubaser'*Xt(:,1));

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
    Ydmd = Ydmd+Phi(1,i)*exp(omegai*tspan)*b(i);
    w = [w imag(omegai)];
end
hold on
plot((1:n-4)*dt,Ydmd(2:end-maxlag-1))

Err = abs(Y{s}(3:end-maxlag)-Ydmd(2:end-maxlag-1))./max(abs(Y{s}(3:end-maxlag)));
AErr = sum(Err)/n



