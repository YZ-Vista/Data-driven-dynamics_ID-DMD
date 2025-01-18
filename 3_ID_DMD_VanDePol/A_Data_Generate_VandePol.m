clc
clear
np = 200; fs = 32; dt = 1/fs;
N = fs*np;
tspan = [0:1/fs:(N-1)/fs];

alpha_R = [0.8 0.9 1.0 1.1 1.2]; wsq_R = [0.8 0.9 1.0 1.1 1.2];

%% Continuous time
s = 1;
for i = 1:length(alpha_R)
    for j = 1:length(wsq_R)
        alpha = alpha_R(i);
        wsq = wsq_R(j);
        x0 = [0.1 0];
        [t,x] = ode45(@odetestV,tspan,x0,[],alpha,wsq);
        y = x(:,1).';
        Yvan{s} = y;
        PARA{s} = [alpha wsq];
        s = s+1;
    end
end

save VandePol_data.mat Yvan PARA

