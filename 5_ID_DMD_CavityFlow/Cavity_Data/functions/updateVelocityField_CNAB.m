function [u,v,p] = updateVelocityField_CNAB(u,v,nx,ny,dx,dy,Re,dt,velbc,method,uForce,vForce,velbc1)
% Copyright 2020 The MathWorks, Inc.

persistent Nu_old Nv_old

if isempty(Nu_old) || isempty(Nv_old) || any(size(Nu_old) ~= [nx-1,ny])
    Nu_old = zeros(nx-1,ny);
    Nv_old = zeros(nx,ny-1);
end

% Apply boundary conditions:
% represented as the values on ghost cells
u(:,1) = 2*velbc.uBottom - u(:,2); v(:,1) = velbc.vBottom;             %bottom
u(:,end) = 2*velbc.uTop - u(:,end-1);  v(:,end) = velbc.vTop;  %top
u(1,:) = velbc.uLeft;  v(1,:) = 2*velbc.vLeft - v(2,:);             %left
u(end,:) = velbc.uRight;  v(end,:) = 2*velbc.vRight - v(end-1,:);    %right

%  �g�U��(u) Get viscous terms for u
Lux = (u(1:end-2,2:end-1)-2*u(2:end-1,2:end-1)+u(3:end,2:end-1))/dx^2; % nx-1 * ny
Luy = (u(2:end-1,1:end-2)-2*u(2:end-1,2:end-1)+u(2:end-1,3:end))/dy^2; % nx-1 * ny

% �g�U��(v) Get viscous terms for v
Lvx = (v(1:end-2,2:end-1)-2*v(2:end-1,2:end-1)+v(3:end,2:end-1))/dx^2; % nx * ny-1
Lvy = (v(2:end-1,1:end-2)-2*v(2:end-1,2:end-1)+v(2:end-1,3:end))/dy^2; % nx * ny-1

% �Η����̌v�Z
% Get nonlinear terms
% 1. interpolate velocity at cell center/cell cornder
uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
uco = (u(:,1:end-1)+u(:,2:end))/2;
vco = (v(1:end-1,:)+v(2:end,:))/2;
vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;

% 2. multiply
uuce = uce.*uce;
uvco = uco.*vco;
vvce = vce.*vce;

% 3-1. get derivative for u
Nu = (uuce(2:end,:) - uuce(1:end-1,:))/dx;
Nu = Nu + (uvco(2:end-1,2:end) - uvco(2:end-1,1:end-1))/dy;
% forcing
Nu = Nu - uForce;

% 3-2. get derivative for v
Nv = (vvce(:,2:end) - vvce(:,1:end-1))/dy;
Nv = Nv + (uvco(2:end,2:end-1) - uvco(1:end-1,2:end-1))/dx;
% forcing
Nv = Nv - vForce;


%  ���̑��x��Z�o
% Implicit treatment for xy direction
Lubc = zeros(size(Luy));
Lubc(1,:) = velbc1.uLeft(2:end-1)/dx^2;
Lubc(end,:) = velbc1.uRight(2:end-1)/dx^2;
Lubc(:,1) = 2*velbc1.uBottom(2:end-1)/dy^2;
Lubc(:,end) = 2*velbc1.uTop(2:end-1)/dy^2;

Lvbc = zeros(size(Lvy));
Lvbc(1,:) = 2*velbc1.vLeft(2:end-1)/dx^2;
Lvbc(end,:) = 2*velbc1.vRight(2:end-1)/dx^2;
Lvbc(:,1) = velbc1.vBottom(2:end-1)/dy^2;
Lvbc(:,end) = velbc1.vTop(2:end-1)/dy^2;

b = u(2:end-1,2:end-1) - dt*((3*Nu-Nu_old)/2 - 1/(2*Re)*(Lux+Luy+Lubc));
xu = getIntermediateU_xyCNAB(u, b, dt, Re, nx, ny, dx, dy);
b = v(2:end-1,2:end-1) - dt*((3*Nv-Nv_old)/2 - 1/(2*Re)*(Lvx+Lvy+Lvbc));
xv = getIntermediateV_xyCNAB(v, b, dt, Re, nx, ny, dx, dy);

u(2:end-1,2:end-1) = xu;
v(2:end-1,2:end-1) = xv;

Nu_old = Nu;
Nv_old = Nv;

% represented as the values on ghost cells
u(:,1) = 2*velbc1.uBottom - u(:,2); v(:,1) = velbc1.vBottom;             %bottom
u(:,end) = 2*velbc1.uTop - u(:,end-1);  v(:,end) = velbc1.vTop;  %top
u(1,:) = velbc1.uLeft;  v(1,:) = 2*velbc1.vLeft - v(2,:);             %left
u(end,:) = velbc1.uRight;  v(end,:) = 2*velbc1.vRight - v(end-1,:);    %right

% �V�������x��
% ���͂̎��i�|���\���������j�������đ��x������ʕۑ��𖞂�����Ɏʑ��B
% velocity correction
% RHS of pressure Poisson eq.
b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
    + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy);

% Solve for p
switch method
    case 'direct'
        % by directly inverting the matrix�i���ږ@�j
        dp = solvePoissonEquation_direct(b,nx,ny,dx,dy);
    case 'dct'
        % by using the discrete cosine transform�i�R�T�C���ϊ��g�p�j
        % Note: Signal Processing Toolbox required
        dp = solvePoissonEquation_2dDCT(b,nx,ny,dx,dy);
    case 'dct_p'
        % by using the discrete cosine transform�i�R�T�C���ϊ��g�p�j
        % Note: Signal Processing Toolbox required
        dp = solvePoissonEquation_2dDCT_p(b,nx,ny,dx,dy);
    
    otherwise
        error("Specified method: " + method + " is not supported." + ...
            "It should be either direct or dct");
end

switch method
    case 'dct_p'
        u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (dp(2:end,:)-dp(1:end-1,:))/dx;
        v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (dp(:,2:end)-dp(:,1:end-1))/dy;
        u(end,2:end-1) = u(end,2:end-1) + 2*dp(end,:)/dx;    %right
        v(2:end-1,end) = v(2:end-1,end) + 2*dp(:,end)/dy;    %top
        v(2:end-1,1) = v(2:end-1,1) - 2*dp(:,1)/dy;    %bottom
    otherwise
        u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (dp(2:end,:)-dp(1:end-1,:))/dx;
        v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (dp(:,2:end)-dp(:,1:end-1))/dy;
end

% correction to get the final velocity
p = dp;

% check the divergence
%     b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
%         + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy);
%     norm(b)

end