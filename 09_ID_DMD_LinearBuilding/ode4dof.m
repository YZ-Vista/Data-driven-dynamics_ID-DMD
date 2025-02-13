function xdot=ode4dof(t,x,tspan,cC1,cC2,cB,cA1,kC1,kC2,kB,kA1)
xdot=zeros(8,1);

mC1=5*10^6; mC2=4*10^6; mA2=3*10^6; mA1=2*10^6;

m1=mC1; m2=mC2; m3=mA2; m4=mA1;
k1=kC1; k2=kC2; k3=kB; k4=kA1;
c1=cC1; c2=cC2; c3=cB; c4=cA1;

M=[m1 0 0 0;
   0 m2 0 0;
   0 0 m3 0;
   0 0 0 m4];
C=[c1+c2 -c2 0 0;
   -c2 c2+c3 -c3 0;
   0 -c3 c3+c4 -c4;
   0 0 -c4 c4];
K=[k1+k2 -k2 0 0;
   -k2 k2+k3 -k3 0;
   0 -k3 k3+k4 -k4;
   0 0 -k4 k4];
xdot(1:4)=x(5:end);
xdot(5:end)=inv(M)*(-C*x(5:end)-K*x(1:4));
end