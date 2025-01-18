clc
clear
np = 200; fs = 32; dt = 1/fs;
N = fs*np;
tspan=[0:1/fs:(N-1)/fs];

k1 = 1e2; c3_set = [1 5 10 15 20];

%% Continuous time
x0=[0.01 0];
for i = 1:length(c3_set)
    c3 = c3_set(i);
[t,x]=ode45(@odetestD,tspan,x0,[],tspan,k1,c3);
y=x(:,1);

Y{i} = y.';
PARA{i} = c3;
end

save Duff_D.mat Y PARA
