function [Phi, Lambda, b, Ubaser, Mcom] = DMD_for_D1(Xsec,Xsec_prime,P,r,pt)
X0=[];X1=[];Xprime=[];
[~,N] = size(Xsec);
for i = 1:N
    X0 = [X0 Xsec{i}];
    X1 = [X1 P{i}*Xsec{i}];
    Xprime = [Xprime Xsec_prime{i}];
end
X=[X0; X1];
[M,~] = size(X0);
% %% SVD 
[Ubase,Sigmabase,Vbase] = svd(Xprime,'econ');
assignin('base', 'Ubase', Ubase'*Ubase);
Ubaser = Ubase(:,1:r);
Sigmabaser = Sigmabase(1:r,1:r);
Vbaser = Vbase(:,1:r);

[U,Sigma,V] = svd(X,'econ');      % Step 1
% assignin('base', 'M', M);
rp = r;
Ur = U(:,1:rp);
UA = U(1:M,:); 
UrA = U(1:M,1:rp);
UB = U(M+1:end,:);
UrB = U(M+1:end,1:rp);
Sigmar = Sigma(1:rp,1:rp);
Vr = V(:,1:rp);

% Atilde = Ubaser'*Xprime*Vr/Sigmar*UrA'*Ubaser;    % Step 2
% Btilde = Ubaser'*Xprime*Vr/Sigmar*UrB'*Ubaser;    % Step 2

Atilde = Ubaser'*Xprime*Vr/Sigmar*UrA'*Ubaser;    % Step 2
Btilde = Ubaser'*Xprime*Vr/Sigmar*UrB'*Ubaser;    % Step 2

CombM = Atilde+pt.*Btilde;

[W,Lambda] = eig(CombM);         % Step 3
% [W,Lambda] = eigs(CombM,120);

%% Projected DMD
% Phi = Ubaser*W;       % Step 4
% b = Phi\X0(:,1);
% Phi =  Xprime*(Vr/Sigmar)*(UrA'+pt.*UrB')*Ubaser*W;       % Step 4
% b = (W*Lambda)\(Ubaser'*X0(:,1));

Phi =  Ubaser*W;       % Step 4
b = (W*Lambda);

%% Reconstuction
% [U,Sigma,V] = svd(X,'econ');
% Acom = Xprime*V/Sigma*U';
% McomC = Ubaser'*(Acom(:,1:M)+pt.*Acom(:,M+1:end))*Ubaser;
% Mcom = Ubaser*McomC*Ubaser';
% Mcom = Acom(:,1:M)+pt.*Acom(:,M+1:end);


Mcom = Ubaser*(Atilde+pt*Btilde)*Ubaser'; 
% 
% % [Phi_M,Lambda_M] = eig(Mcom);
% % b_M = Phi_M\X0(:,1);
% [Phi,Lambda] = eig(Mcom);
% b = Phi\X0(:,1);