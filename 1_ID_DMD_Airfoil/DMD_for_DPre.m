function [Atilde,Btilde,Ubaser] = DMD_for_DPre(Xsec,Xsec_prime,P,r)
X0=[];X1=[];Xprime=[];
[~,N] = size(Xsec);
[~,Nt] = size(Xsec{1});
for i = 1:N
    X0 = [X0 Xsec{i}];
    X1 = [X1 P{i}*Xsec{i}];
    Xprime = [Xprime Xsec_prime{i}];
end
X = [X0; X1];
[M,~] = size(X0);

%     randi('state',tr)
bag_idx = sort(randperm(N*Nt,round(N*Nt/2)));
X_bag = X(:,bag_idx);
Xprime_bag = Xprime(:,bag_idx);
% %% SVD
[Ubase,Sigmabase,Vbase] = svd(Xprime_bag,'econ');
assignin('base', 'Ubase', Ubase'*Ubase);
Ubaser = Ubase(:,1:r);
Sigmabaser = Sigmabase(1:r,1:r);
Vbaser = Vbase(:,1:r);

[U,Sigma,V] = svd(X_bag,'econ');      % Step 1
% assignin('base', 'M', M);
rp = r;
Ur = U(:,1:rp);
UA = U(1:M,:);
UrA = U(1:M,1:rp);
UB = U(M+1:end,:);
UrB = U(M+1:end,1:rp);
Sigmar = Sigma(1:rp,1:rp);
Vr = V(:,1:rp);

Atilde = Ubaser'*Xprime_bag*Vr/Sigmar*UrA'*Ubaser;    % Step 2
Btilde = Ubaser'*Xprime_bag*Vr/Sigmar*UrB'*Ubaser;    % Step 2
