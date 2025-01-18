clc
clear
load BUILDING;
X = BUILDING;
r = 16 %
indx=[1 3 5 6];
for i = 1:length(indx)
    k=indx(i);
    Xc{i} = X{k}(:,1:end-1);
    Xc_prime{i} = X{k}(:,2:end);
    P{i} = PARA{k};
end
s = 1;
[m,n] = size(X{s});

for run = 1:60
    pt_min = 2;
    pt_max = 3;
    fs = 80; dt=1/fs; np=100; N=np*fs;
    
    for sec = 1:10
        pt = (pt_min+pt_max)/2;
        [Phi, Lambda, b] = DMD_for_D1_bag(Xc,Xc_prime,P,r,pt);
        
        %% Mode prediction
        omega = log(Lambda)/dt;
        tspan = [0:1/fs:(N-1)/fs];
        Y = zeros(1,N);
        w = [];
        for i=1:r
            if real(omega(i,i))>0
                omegai = omega(i,i)-real(omega(i,i));
            else
                omegai = omega(i,i);
            end
            Y = Y+Phi(1,i)*exp(omegai*tspan)*b(i);
            w = [w abs(omegai)];
        end
        Pome = min(w);
        if Pome>=10
            pt_max = pt;
        else
            pt_min = pt;
        end
    end
    Pole(run) = pt;
    idx = run
end

figure
boxplot(Pole);

Spole = std(Pole,1)
Mpole = mean(Pole)


