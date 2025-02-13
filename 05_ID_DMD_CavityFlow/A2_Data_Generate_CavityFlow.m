ReD = [400 500 600 700 800]; % Reynolds number
VetopD = [0.4 0.5 0.6 0.7 0.8]; % top velocity
nt = 3000; % max time steps
Lx = 1; Ly = 1; % domain size
Nx = 80; Ny = 80; % Number of grids
dt = 0.01; % time step;
% Grid size (Equispaced)
dx = Lx/Nx;
dy = Ly/Ny;
% Coordinate of each grid (cell center)
xce = ((1:Nx)-0.5)*dx;
yce = ((1:Ny)-0.5)*dy;
% Coordinate of each grid (cell corner)
xco = (0:Nx)*dx;
yco = (0:Ny)*dy;

s = 1;
for i = 1:length(ReD)
    for j = 1:length(VetopD)
        u = zeros(Nx+1,Ny+2); % velocity in x direction (u)
        v = zeros(Nx+2,Ny+1); % velocity in y direction (v)
        uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center
        vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center

        recordGIF = false; % set to true if you wanna make GIF
        recordRate = 20;
        filename = 'animation_sample.gif'; % Specify the output file name

        figure
        [Xce,Yce] = meshgrid(xce,yce); % cell center grid
        [~,h_abs] = contourf(Xce',Yce',sqrt(uce.^2+vce.^2),'LineColor','none'); % contour

        xlim([0 Lx]); ylim([0 Ly]);

        haxes = gca;
        haxes.XTick = [];
        haxes.YTick = [];

        initialFrame = true;

        Re = ReD(i);
        Ve = VetopD(j); % top velocity
        Tstep = 0;
        for ii = 1:nt
            % Update the velocity field (uses dct)
            [u,v] = updateVelocityField_Euler_bctop(u,v,Nx,Ny,dx,dy,Re,dt,Ve,'dct');
    
            % Update the plot at every recordRate steps
            if mod(ii,recordRate) == 0
                % get velocity at the cell center (for visualization)
                uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center
                vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center

                h_abs.ZData = sqrt(uce.^2+vce.^2)
                h_abs.LevelList = [0:0.01:0.5];
                drawnow
                Tstep = Tstep + 1;
                yout = h_abs.ZData;
                Y_Cavity{s}(:,Tstep) = reshape(yout,Nx*Ny,1);    
            end
        end
        PARA{s} = [Re/1e3 Ve];
        s = s+1;
    end
end
save Cavity_D.mat Y_Cavity PARA Nx Ny