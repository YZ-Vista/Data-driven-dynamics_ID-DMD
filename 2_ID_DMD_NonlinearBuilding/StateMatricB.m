function [X, X_prime] = StateMatricB(Y_data,Y_sec,ny,term_degree)
maxlag = ny;
Nk=size(Y_data,1)-maxlag;  
[~,nm1]=size(Y_data);
%% Construction of regression matrices using data:NARX
no_terms=zeros(1,term_degree); % no of the monomials with each specified degree
no_terms(1)=1; %  no of the 'process' monomials of degree 1
for L=2:term_degree+1
    j=L-1; % degress of the 'process' monomials
    no_terms(L)=no_terms(L-1)*(2*ny+j-1)/j;
end
M=sum(no_terms)-1; 
% total number of candidate terms in the polynomial NARX model without constant
T=zeros(Nk,M+2*ny); % regression matrix (for NARX process model)
T_sec=zeros(Nk,M); % regression matrix (for NARX process model)
for k1 = 1:nm1
    for j = 1:ny
        T(:,ny*(k1-1)+j) = Y_data(maxlag+1-j:maxlag+Nk-j,k1); 
% the monomials of degress 1 (formed by y(k-1)...y(k-ny)): [y(k-1)...y(k-ny)]
    end
end

[~,nm2]=size(Y_sec);
for k2 = 1:nm2
    for j = 1:ny
        T_sec(:,ny*(k2-1)+j) = Y_sec(maxlag+1-j:maxlag+Nk-j,k2); 
% the monomials of degress 1 (formed by y(k-1)...y(k-ny)): [y(k-1)...y(k-ny)]
    end
end

for L=2:term_degree
  jnot=2*ny+sum(no_terms(1:L))-1;
  T(:,jnot+1:jnot+no_terms(L+1)) = State_Polyconstruct(T_sec(:,1:2*ny),L);...
        % the monomial of degree L
end
X = T(1:end-1,:).';
X_prime = T(2:end,:).';
end
