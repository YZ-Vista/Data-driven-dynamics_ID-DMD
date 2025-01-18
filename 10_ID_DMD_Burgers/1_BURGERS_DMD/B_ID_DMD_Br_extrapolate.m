clear all, 
clc
load Map_rebu.mat Map_rebu;
load Map_grye.mat Map_grye;
load Burgers.mat Bur PARA
% VORTALL contains flow fields reshaped into column vectors
X = Bur;
r = 40;
[~,n] = size(X);
XL = X;

[Mx,Nx] = size(X{1}(:,1:101))
EX=[];
A = 0; % no noise
for ni=1:Nx
randn('state',ni)
EX(:,ni) = 2*A*randn(1,Mx)-A;
end
for xi=1:n
    Xdn{xi} = X{xi}(:,1:101);
    X{xi} = Xdn{xi}+EX;
end

tic
indx=[3:4:19];
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:end-1);
    Xc_prime{i} = X{k}(:,2:end);
    P{i} = PARA{k}*1;
    X_train_save{i} = X{k}
end
s = 1; % for pt = 0.01
[Ms,Ns] = size(X{s}(:,1:end-1));
pt = PARA{s}*1;

[Phi, Lambda, b, Ubaser, Mcom] = DMD_for_D1(Xc,Xc_prime,P,r,pt);
% Code to plot the second mode
b = b\(Ubaser'*Xdn{s}(:,1));
toc

Y = zeros(Mx, Nx);
Y(:,1) = Xdn{s}(:,1);

dt=0.01;
tspan = [0:dt:(Nx+200-1)*dt];

Err = zeros(Mx, Nx);
%% Optimised
% Precompute omega and omegai outside the loop
omega = log(Lambda) / dt;
omegai = diag(omega);  % Extract diagonal elements of omega
% Modify omegai based on the real part condition
positive_real_idx = real(omegai) > 0;  % Logical index of elements with positive real part
omegai(positive_real_idx) = omegai(positive_real_idx) - real(omegai(positive_real_idx));

tic
% Compute the exponentials only once for all tspan and omegai
exponent_factors = exp(omegai*tspan); % Precomputes all i and k in one step
% Element-wise scaling b to match the dimensions of `exponent_factors`
scaled_b = exponent_factors.* b;
% Use matrix multiplication for the entire Y computation
Y = Phi*scaled_b;
toc

Y(:,1) = Xdn{s}(:,1);
for k = 1:Nx+200
    Err(:,k) = abs(real(Y(:,k))-XL{s}(:,k))./max(XL{s}(:,k));
end

t = 0:0.01:3;
tT = 0:0.01:1;
x = 0:0.01:1;

figure (1)
subplot(1,3,1)
surf(t(1:3:end),x,real(Y(:,1:3:end)))
shading interp
colormap(subplot(1,3,1),Map_rebu);
caxis([-0.5 0.5])
view([0,90])

subplot(1,3,2)
surf(t(1:3:end),x,XL{s}(:,1:3:end)), 
shading interp
colormap(subplot(1,3,2),Map_rebu);
caxis([-0.5 0.5])
view([0,90])

subplot(1,3,3)
surf(t,x,Err)
shading interp
colormap(subplot(1,3,3),Map_grye);
% caxis([0 6])
view([0,90])

pt = PARA{s};
Burg_pre = real(Y);
Burg_true = XL{s};
% save Bur_Test_Extra.mat x t Burg_pre Burg_true pt