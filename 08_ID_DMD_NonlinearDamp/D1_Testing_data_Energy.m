clc
clear
np = 200; fs = 32; dt = 1/fs;
N = fs*np;
tspan=[0:1/fs:(N-1)/fs];

k1 = 1e2; c3_set = [1:0.1:20];

%% Continuous time
x0=[0.01 0];
c3=0;
[t,x]=ode45(@A_odetestD,tspan,x0,[],tspan,k1,c3);
y0=x(:,1);
Ener0 = sum(y0.^2);
for i = 1:length(c3_set)
    c3 = c3_set(i);
[t,x]=ode45(@odetestD,tspan,x0,[],tspan,k1,c3);
y=x(:,1);

YE{i} = y.';
Ener(i) = (Ener0-sum(y.^2))/Ener0;
PARAE{i} = c3;
end
plot(c3_set,Ener,'.')
hold on

save Duff_D_Energy.mat YE PARAE
