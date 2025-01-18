function T = State_Polyconstruct(X,L)
% This function computes the regression matrix constructed by all the monomials
% of degree L in a multivariate polynomial using the measured data (to be 
% used with Matlab functions 'polyNARMAX', 'polyNARMAXaOLS', 'polyNARMAXpressOLS' 
% developed by Dr. Ping Li)
%
% Inputs of Function:
% X      N by n data matrix [x1 x2 ... xn] where xi (i=1,2,...,n) is a
%        column vector of length N containing the measured value of ith
%        variable (xi).
% L      degree of monomial
%
% Outputs of Function:
% P       N by n_L(=n_{L-1}(n+L-1)/L) regression matrix constructed by the  
%         monomials of degree L in a multivariate polynomial
% p_info  one dimensional cell array of length n_L containing information 
%         on how each of the monomials correponding to a column of
%         the matrix P are made from
%
% Dr Ping Li September 2010
% copyright (c) 2010 Dr Ping Li
[N,n]=size(X);
no_terms=zeros(1,L+1); % no of terms with each of the specified degree
no_terms(1)=1; %  no of terms in a polynomial of degree 0
for J=2:L+1
    j=J-1; % degress of polynomial
    no_terms(J)=no_terms(J-1)*(n+j-1)/j;
end
T=ones(N,no_terms(L+1)); redo=0;
M=sum(no_terms(1:L))-1; % total number of terms in a multivariate polynomial
                          % up to the degress of "L-1"
nL=1; % column No. in P matrix
iL=ones(1,L); % array of length L containing the series numbers of the variable
              % in a monomial of degree L. iL=[i_1=1:n, i_2=i_1:n,..., i_L=i_{L-1}:n] 
while nL<=no_terms(L+1)
    for j=1:L
        if iL(j)<=n
            T(:,nL)=T(:,nL).*X(:,iL(j));
            if j==L
%                 p_info{nL}=struct('column_No_in', nL+M, 'variable_No_ny_nu_np', iL);
                iL(L)=iL(L)+1;
            end
        elseif j>1
            iL(j-1)=iL(j-1)+1;
            for ij=j:L
                iL(ij)=iL(ij-1);
            end
            T(:,nL)=ones(N,1); redo=1; break
        end
    end
    if redo>0
        redo=0;
    else
        nL=nL+1;
    end
end
% end of the function