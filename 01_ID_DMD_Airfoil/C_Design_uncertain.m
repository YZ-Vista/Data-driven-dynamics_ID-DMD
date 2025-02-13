clear all, 
clc
%% Prediction
load Ene.mat Eness EnessRate Theta

up = max(Eness)-mean(Eness);
down = mean(Eness)-min(Eness);
errorbar(Theta,mean(Eness),down,up,'.-');
hold on

load VORTEX.mat VORTEX PARA CYLIN lx ly;
PARA_D = PARA;
[M,N] = size(VORTEX{1});
EX=[];
A = 0e-1; % 3e-1
for ni=1:N
randn('state',ni)
EX(:,ni) = 2*A*randn(1,ly*lx)-A;
end
for sa = 1:11
  X{sa} = VORTEX{sa}(:,40:end)+EX(:,40:end);
end
%%
% VORTALL contains flow fields reshaped into column vectors
vote = [1:11];
for cd = 1:length(vote)
s = vote(cd);
Theta(s) = PARA_D{s}*180/pi;
[M,N] = size(X{s});

for k=1:N
    Mflow(:,:,k) = reshape(X{s}(:,k),ly,lx);
end
[L,W,T] = size(Mflow);

for i = 1:N
EneM(i) = sum(sum(Mflow(:,145:end,i).^2));
end
SEness(s) = abs(sqrt(mean(EneM)));
end

plot(Theta,SEness,'o')
hold on
