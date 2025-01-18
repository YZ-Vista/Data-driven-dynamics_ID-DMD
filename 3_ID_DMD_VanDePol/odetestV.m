function xdot=odetestV(t,x,alpha,wsq)
% Van de Pol equation: x'' -alpha(1-x^2)*x'+w^2*x=0
xdot=zeros(2,1);
xdot(1)=x(2);
xdot(2)=alpha*(1-x(1)^2)*x(2)-wsq*x(1);
end

