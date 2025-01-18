clc
clear
load VandePol_data.mat Yvan PARA;

[~,nm] = size(Yvan);
for xi = 1:nm
    XS{xi} = Yvan{xi}(:,1000:end).';
    ntrunc = round(1*length(XS{1}));
    X{xi} = XS{xi}(1:ntrunc,:);
end
N = round(length(XS{1}));

indx = [1:2:25]; %
term_degree = 11;
ny = 2; % Time delay for output y (k-1,k-2);
maxlag = ny;
for i = 1:length(indx)
    k=indx(i);
    [Xc{i}, Xc_prime{i}] = StateMatric(X{k},ny,term_degree);
    P{i} = PARA{k};
end
s = 14; % 2
[Xt, ~] = StateMatric(X{s},ny,term_degree);
pt = PARA{s};
%% DMD for design
r = 70;
[Phi, Lambda, b, Ubaser] = DMD_for_D1(Xc,Xc_prime,P,r,pt);

fs = 32; dt = 1/fs;
plot(XS{s}(maxlag+1:end-2),(XS{s}(maxlag+1:end-2)-XS{s}(maxlag:end-3))/dt)
hold on

XPm(:,1)=Xt(:,1);
b=b\(Ubaser'*Xt(:,1));

%% Mode prediction
omega = log(Lambda)/dt;
tspan = [0:1/fs:(N-1)/fs];
Ydmd = zeros(1,N);
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

plot(Ydmd(2:end-maxlag-1),(Ydmd(2:end-maxlag-1)-Ydmd(1:end-maxlag-2))/dt)

Err = abs(XS{s}(maxlag+1:end-2).'-Ydmd(2:end-maxlag-1))./max(abs(XS{s}(maxlag+1:end-2)));

AErr = sum(Err)/N


