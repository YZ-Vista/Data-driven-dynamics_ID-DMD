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
%% Testing data
load Duff_D_Energy.mat YE PARAE;
fs = 32; dt = 1/fs;
tspan=[0:1/fs:(n-1)/fs];

k1 = 1e2;c3=0;
x0=[0.01 0];

[t,x]=ode45(@odetestD,tspan,x0,[],tspan,k1,c3);
y0=x(:,1);
Ener0 = sum(y0.^2);

for s = 1:length(PARAE); % 2
    [Xt, ~] = StateMatric(YE{s},ny,term_degree);
    pt(s) = PARAE{s}/10;
    
    [t,x]=ode45(@odetestD,tspan,x0,[],tspan,k1,PARAE{s});
    ysim=x(:,1);
    EnerSim(s) = real((Ener0-sum(ysim.^2))/Ener0);
    %% DMD for design
    r = 30;
    %% Find Uncertainty
    for run = 1:100
        [Phi, Lambda, b, Ubaser] = DMD_for_D1_bag(Xc,Xc_prime,P,r,pt(s));
        
        %% Plot time series
        XPm(:,1)=Xt(:,1);
        b=b\(Ubaser'*Xt(:,1));
        
        %% Mode prediction
        omega = log(Lambda)/dt;
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
        
        Ener(run,s) = real((Ener0-sum(Ydmd.^2))/Ener0);
    end
    %% Design results
    [Phi, Lambda, b, Ubaser] = DMD_for_D1(Xc,Xc_prime,P,r,pt(s));

    %% Plot time series
    XPm(:,1)=Xt(:,1);
    b=b\(Ubaser'*Xt(:,1));

    %% Mode prediction
    omega = log(Lambda)/dt;
    Ydmdm = zeros(1,n);
    w = [];
    for i=1:r
        if real(omega(i,i))>0
            omegai = omega(i,i)-real(omega(i,i));
        else
            omegai = omega(i,i);
        end
        Ydmdm = Ydmdm+Phi(1,i)*exp(omegai*tspan)*b(i);
        w = [w imag(omegai)];
    end
    Enerm(s) = real((Ener0-sum(Ydmdm.^2))/Ener0);
    s
end
% save Ener.mat Ener pt
Ener = rmoutliers(Ener,1);
xconf = [pt pt(end:-1:1)] ;% forward-backward of x-axis
Eup = max(Ener,[],1);
Elow = min(Ener,[],1);
yconf = [Eup Elow(end:-1:1)];% upper and lower bound

figure
p = fill(xconf,yconf,'red');% uncertain interval
p.FaceColor = [0.6 0.8 0.8];% colour    
p.EdgeColor = 'none';% Bound colour

hold on
plot(pt,EnerSim,'r')
hold on
plot(pt,Enerm,'b')
