function [Theta, Eness] = Uncertain_Design(VORTEX_D, PARA_D, lx, ly)
X = VORTEX_D;
% VORTALL contains flow fields reshaped into column vectors
vote = [1:11];
for cd = 1:length(vote)
s = vote(cd);
Theta(s) = PARA_D{s}*180/pi;
[~,N] = size(X{s});

dt=1/36;
tspan = [dt:dt:N*dt];

for k=1:N
    Mflow(:,:,k) = reshape(X{s}(:,k),ly,lx);
end

for i = 1:N
EneM(i) = sum(sum(Mflow(:,145:end,i).^2));
end
Eness(s) = abs(sqrt(mean(EneM)));
end

% plot(Theta,Eness, 'Color', [0.7, 0.7, 0.7])
% hold on