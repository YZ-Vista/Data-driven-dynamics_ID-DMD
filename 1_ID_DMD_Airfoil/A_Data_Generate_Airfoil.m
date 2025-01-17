clc
clear
%% Note: It will take 24h to generate all the data via simulation
load Flow.mat Flow;
Nx = 120*6; lx = Nx+2; 
Ny =  50*6; ly = Ny+2;
Re = 150;
Lx =   12;
Ly =   5;
dt = 0.005;
Ulid = 1;

dx = Lx/Nx;
dy = Ly/Ny;
Theta_set = -[2:1:12]*pi/180; % Different pitched angle
tic
for s = 1:length(Theta_set)
u = ones(Ny+2, Nx+2); v= zeros(Ny+2, Nx+2); 
p = u; Dp = p;
U = ones(Ny, Nx+1);
V = zeros(Ny+1, Nx);
Nmax = 10000;
dts = zeros(Nmax,1);
Fs  = zeros(Nmax,2);
% cfl = 0.2;
ipx = 2:(Nx+1); % collocated; phyiscal
ipy = 2:(Ny+1); % collocated; physical
us = u; % u star; physical+ghost
vs = v; % v star; physical+ghost
Us = U; Vs = V;   % Us and Vs for boundaries

% setup for poisson
kx = (0 : Nx-1);
ky = (0 : Ny-1);
[kx,ky] = meshgrid(kx,ky);
mwx = 2 * (cos(pi * kx / Nx)-1)/dx^2;
mwy = 2 * (cos(pi * ky / Ny)-1)/dy^2;

%% For the cylinder
xx = linspace(-dx/2, Lx+dx/2,Nx+2);
yy = linspace(-dy/2, Ly+dy/2,Ny+2);
[x, y] = meshgrid(xx, yy);
cx = Lx/8;
cy = Ly/2;
r = 0.3;
ra = 0.6;
rb = 0.3;

theta = Theta_set(s);
% Define the shape of the airfoil
isInCircle = ((x-cx)*cos(theta)-(y-cy)*sin(theta)).^2 /ra^2 ...
    + ((x-cx)*sin(theta)-(y-cy)*cos(theta)).^2 /rb^2 <= 1;  % same size as u and v
onBounary = isInCircle;
onBounary(ipy,ipx) = isInCircle(ipy,ipx) & (...
                      xor(isInCircle(ipy,ipx+1), isInCircle(ipy,ipx-1)) |...
                      xor(isInCircle(ipy+1,ipx), isInCircle(ipy-1,ipx))   );   
boundL_neigh = false(Ny+2, Nx+2);
boundR_neigh = boundL_neigh; boundT_neigh = boundL_neigh; boundB_neigh = boundL_neigh;
boundL_neigh(ipy,ipx) =  onBounary(ipy,ipx+1);
boundR_neigh(ipy,ipx) =  onBounary(ipy,ipx-1);
boundT_neigh(ipy,ipx) =  onBounary(ipy-1,ipx);
boundB_neigh(ipy,ipx) =  onBounary(ipy+1,ipx);
InnerCircle = isInCircle & (~ onBounary);

%% for implicit solver
mu = 1/Re;
e = ones(Nx,1);
Ax = spdiags([-mu*0.5*dt*e./dx^2 (1+mu*dt/dx^2)*e -mu*0.5*dt*e./dx^2],-1:1,Nx,Nx); 
e = ones(Ny,1);
Ay = spdiags([-mu*0.5*dt*e./dy^2 (1+mu*dt/dy^2)*e -mu*0.5*dt*e./dy^2],-1:1,Ny,Ny); 

k = 0;
for i = 1:Nmax
% interplation from u->U
u_w = (u(ipy, ipx-1)+u(ipy, ipx))/2; u_e = (u(ipy, ipx+1)+u(ipy, ipx))/2;
u_n = (u(ipy+1, ipx)+u(ipy, ipx))/2; u_s = (u(ipy-1, ipx)+u(ipy, ipx))/2;

v_n = (v(ipy+1, ipx)+v(ipy, ipx))/2; v_s = (v(ipy-1, ipx)+v(ipy, ipx))/2;
v_w = (v(ipy, ipx-1)+v(ipy, ipx))/2; v_e = (v(ipy, ipx+1)+v(ipy, ipx))/2;

conu = (U(:, 2:Nx+1).*u_e - U(:, 1:Nx).*u_w)/dx + ...
       (V(2:Ny+1,: ).*u_n - V(1:Ny, :).*u_s)/dy;

conv = (U(:, 2:Nx+1).*v_e - U(:, 1:Nx).*v_w)/dx + ...
       (V(2:Ny+1,: ).*v_n - V(1:Ny, :).*v_s)/dy;

us(ipy, ipx) = u(ipy,ipx) - dt * conu;
vs(ipy, ipx) = v(ipy,ipx) - dt * conv;
us = ADI(dt, dx, dy, 1/Re, us, Ax, Ay, ipy, ipx,Ny, Nx);
vs = ADI(dt, dx, dy, 1/Re, vs, Ax, Ay, ipy, ipx,Ny, Nx);

% interpolation
Us(:, :) = (us(ipy,1:Nx+1)+us(ipy,2:Nx+2))/2;
Vs(:, :) = (vs(1:Ny+1,ipx)+vs(2:Ny+2,ipx))/2;

Us(1:Ny, 1) =    Ulid; %left bc
Us(1:Ny, Nx+1) = Ulid; %Right bc

F = ((Us(:,2:Nx+1) - Us(:,1:Nx))/dx + (Vs(2:Ny+1,:) - Vs(1:Ny,:))/dy)/dt;
F(isInCircle(ipy,ipx)) = 0;

%% Solve the poisson equation
B = dct2(F);
A =  B ./ ( mwx+mwy ); A(1,1) = 0;
p(ipy, ipx) = idct2(A);

% BC for p:
p(:, 1) = p(:,2);
p(:,Nx+2) = p(:, Nx+1);
p(1,:) = p(2,:);
p(Ny+2, :) = p(Ny+1,:);

%% The third step:
% for cell centered values:
u(ipy, ipx) = us(ipy, ipx) - dt * (p(ipy,ipx+1)-p(ipy,ipx-1))/(2*dx);
v(ipy, ipx) = vs(ipy, ipx) - dt * (p(ipy+1,ipx)-p(ipy-1,ipx))/(2*dy);
% for face centered values:
U(:, 2:Nx) =   Us(:, 2:Nx) -dt * (p(ipy, 3:Nx+1) - p(ipy, 2:Nx))/dx;
V(2:Ny, :) =   Vs(2:Ny, :) -dt * (p(3:Ny+1,ipx) -  p(2:Ny,ipx))/dy;

% Impose boundary for u:
u(1,:) =       u(2,:);               %Down D
u(:, Nx+2) =  2*Ulid - u(:, Nx+1);   %Right R
u(:, 1) =     2*Ulid -  u(:, 2);     %Left L
u(Ny+2, :) =     u(Ny+1,:);          %Up U

% Impose boundary for v:
v(:, 1) = -v(:, 1);                   %Left L
v(:,Nx+2) =  v(:,Nx+1);               %Right R
v(1,:) =  -v(2,:);                    %Down D
v(Ny+2, :) = -v(Ny+1, :);             %Up U

% % for cylinder:
u(isInCircle) = 0; v(isInCircle) = 0;

Fs(i, 1) =  2*dx*sum( 0.5*(p(onBounary)+p(boundL_neigh)).*isInCircle(boundR_neigh) +...
                   -  0.5*(p(onBounary)+p(boundR_neigh)).*isInCircle(boundL_neigh));
Fs(i, 2) = -2*dy*sum( 0.5*(p(onBounary)+p(boundT_neigh)).*isInCircle(boundB_neigh) +...
                  -   0.5*(p(onBounary)+p(boundB_neigh)).*isInCircle(boundT_neigh));

if mod(i,100)==0
    [dux, duy] = gradient(u);
    [dvx, dvy] = gradient(v);
    vort = duy/(2*dy) - dvx/(2*dx);
    vort(isInCircle) = 0;
% Plot
    contourf(vort,-20:0.05:20, 'LineColor','none')
    title(['t=',num2str(dt*i)])
    colormap(Flow);
    caxis([-2 2]) % clim after MATLAB 2022
    colorbar;
    axis equal
    drawnow;
    k = k+1;
    VOE(:,k) = reshape(vort,lx*ly,1); 
end
end
    VORTEX{s} = VOE;
    PARA{s} = -theta;
    CYLIN{s} = isInCircle;
end
save VORTEX.mat VORTEX PARA CYLIN lx ly
toc
div = (v(ipy+1, ipx)-v(ipy-1, ipx))/(2*dy)+...
      (u(ipy, ipx+1)-u(ipy, ipx-1))/(2*dx);


function u = ADI(dt, dx, dy, mu, u, Ax, Ay, ipy, ipx,Ny,Nx)
% %------------implicit in x direction-----------
b =  u(ipy,ipx)+0.5 * dt * mu * (u(ipy+1, ipx)+u(ipy-1,ipx)-2*u(ipy,ipx))/dy^2;
b(ipy-1,1  ) = b(ipy-1,1 )  +mu*0.5*dt*u(ipy,    1)/dx^2;
b(ipy-1,Nx ) = b(ipy-1,Nx)  +mu*0.5*dt*u(ipy, Nx+2)/dx^2;
u(ipy,ipx) = (Ax\(b(ipy-1, :).'))';

% % -----------implicit in y direction--------
b =  u(ipy,ipx)+0.5 * dt * mu * (u(ipy, ipx+1)+u(ipy,ipx-1)-2*u(ipy,ipx))/dx^2;
b(1,ipx-1 )= b(1  ,ipx-1)+   mu*0.5*dt*u(1,    ipx)/dy^2;
b(Ny,ipx-1)= b(Ny,ipx-1) +   mu*0.5*dt*u(Ny+2, ipx)/dy^2;
u(ipy,ipx) = Ay\b(:,ipx-1);
end
