function xdot=A_odetestT(t,x,omega,k1,c3)
% Duffing equation: x'' +0.03x'+1e2*x+c3*x'^3=0
u=cos(omega*t);
xdot=zeros(2,1);
xdot(1)=x(2);
xdot(2)=u-k1*x(1)-0.1*x(2)-c3*x(2)^3;
end

