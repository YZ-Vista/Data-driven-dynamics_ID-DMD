clc
clear
addpath('D:\Backup\Working papers\2024 DMD for design\GitHub_upload\10_ID_DMD_Burgers\1_BURGERS\chebfun-chebfun-2b127f1'), savepath
% number of realizations to generate
N = 1;

% parameters for the Gaussian random field
gamma = 4;
tau = 5;
sigma = 25^(2);

% viscosity
visc_set = [0.01:0.002:0.05];

for z = 1:length(visc_set)
    visc = visc_set(z);
% grid size
    s = 4096;
    steps = 300;
    nn = 101;

    input = zeros(N, nn);
    if steps == 1
        output = zeros(N, s);
    else
        output = zeros(N, steps, nn);
    end

    tspan = linspace(0,3,steps+1);
    x = linspace(0,1, s+1);
    X = linspace(0,1, nn);

    for j=1:N
        u0 = GRF(s/2, 0, gamma, tau, sigma, "periodic");
        u = Burgers(u0, tspan, s, visc);

        u0_eval = u0(X);
        input(j,:) = u0_eval;

        if steps == 1
            output(j,:) = u.values;
        else
            for k=1:(steps+1)
                output(j,k,:) = u{k}(X);
            end
        end

        disp(j);
    end

    Bur{z} = squeeze(output(1,:,:)).';
    PARA{z} = visc;
    BurT{z} = output(1,:,:);
end
save Burgers.mat Bur PARA
t = 0:0.01:3;
x = 0:0.01:1; 
surf(t,x,Bur{3})
view([0,90])
