function [Phi, Lambda, b] = DMD_for_D1_bag(Xsec,Xsec_prime,P,r,pt)
X0=[];X1=[];Xprime=[];
[~,N] = size(Xsec);
[~,Nt] = size(Xsec{1});
for i = 1:N
    X0 = [X0 Xsec{i}];
    X1 = [X1 P{i}*Xsec{i}];
    Xprime = [Xprime Xsec_prime{i}];
end
X=[X0; X1];
[M,~] = size(X0);

bag_idx = sort(randperm(N*Nt,round(N*Nt/2)));
X_bag = X(:,bag_idx);
Xprime_bag = Xprime(:,bag_idx);

%% SVD 
[Ubase,Sigmabase,Vbase] = svd(Xprime_bag,'econ');
Ubaser = Ubase(:,1:r);
Sigmabaser = Sigmabase(1:r,1:r);
Vbaser = Vbase(:,1:r);

[U,Sigma,V] = svd(X_bag,'econ');      % Step 1
rp = r;
Ur = U(:,1:rp);
UrA = U(1:M,1:rp);
UrB = U(M+1:end,1:rp);
Sigmar = Sigma(1:rp,1:rp);
Vr = V(:,1:rp);

Atilde = Ubaser'*Xprime_bag*Vr/Sigmar*UrA'*Ubaser;    % Step 2
Btilde = Ubaser'*Xprime_bag*Vr/Sigmar*UrB'*Ubaser;    % Step 2

CombM = Atilde+pt.*Btilde;

[W,Lambda] = eig(CombM);         % Step 3

%% Projected DMD
% Phi = Ubaser*W;       % Step 4
% b = Phi\X0(:,1);
Phi =  Xprime_bag*(Vr/Sigmar)*(UrA'+pt.*UrB')*Ubaser*W;       % Step 4
b = (W*Lambda)\(Ubaser'*X0(:,1));
