clc
clear
% Parameters
Lx = 5.0;    % x length
Ly = 5.0;    % y length
Nx = 50;     % x mesh number
Ny = 50;     % y mesh number
dx = Lx / (Nx - 1); % x mesh step
dy = Ly / (Ny - 1); % y mesh step
alpha_set = [0.2:0.1:0.8]; % Thermal conductivity
dt = 0.0005;  % Time step
Nt = 4000;     % Time numbers

for s = 1:length(alpha_set)
    alpha = alpha_set(s);
% Initial thermal matrix
u = zeros(Nx, Ny);
u_new = u;

% Time progress
Tstep = 0;
for n = 1:Nt
    for Ii = 1:3
        for Ij = 1:3
          u(Ii,Ij) = 1;%sin(1*(n-1)*dt);
        end
    end
    % Finite difference method
    for i = 2:Nx-1
        for j = 2:Ny-1
            u_new(i,j) = u(i,j) - 3 * dt * ( ...
                (u(i+1,j) - u(i,j)) / dx + ...
                (u(i,j+1) - u(i,j)) / dy )+ ...
                 alpha * dt * ( ...
                (u(i+1,j) - 2*u(i,j) + u(i-1,j)) / dx^2 + ...
                (u(i,j+1) - 2*u(i,j) + u(i,j-1)) / dy^2 );
        end
    end
    % Boundary condition (i.e. Boundary temperature 0)
    u_new(1,:) = 0;
    u_new(Nx,:) = 0;
    u_new(:,1) = 0;
    u_new(:,Ny) = 0;
    
    % temperature distribution
    u = u_new;
    
    % Plot temperature distribution
    if mod(n, 50) == 0
        Tstep = Tstep + 1;
        Y_Jet{s}(:,Tstep) = reshape(u,Nx*Ny,1); 
        
        surf(linspace(0, Lx, Nx), linspace(0, Ly, Ny), u);
        title(['Time step: ', num2str(n)]);
        xlabel('x');
        ylabel('y');
        zlabel('Temperature');
        view([90,-90])
        drawnow;
    end
    PARA{s} = [alpha];
end
end

save Jet_D.mat Y_Jet PARA Nx Ny
