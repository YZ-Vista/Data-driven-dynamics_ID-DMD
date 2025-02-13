clc
clear
load Duff_D.mat Y PARA;

X = Y;
N = length(X{1});

indx = [1 3 5]; %
term_degree = 8;
ny = 2; % Time delay for output y (k-1,k-2);
maxlag = ny;
for i = 1:length(indx)
    k=indx(i);
    [Xc{i}, Xc_prime{i}] = StateMatric(X{k},ny,term_degree);
    P{i} = PARA{k}/10;
end

pt_T = [1:0.5:20]/10;
for i = 1:length(pt_T)
%% DMD for design
pt = pt_T(i);
[m,n] = size(Xc{1});
r = 35;
    for run = 1:1 % we can consider multiple runs here for bagging
    [Phi, Lambda, b, Ubaser] = DMD_for_D1_bag(Xc,Xc_prime,P,r,pt);

    %% Code to plot the second mode 
    fs = 32; dt = 1/fs;
    [Olam,indf] = sort(diag(log(Lambda)/dt));
    s = 5; % 5, 9, 13
    sec = indf(s);
    % figure
    Vmode = [0; abs(Phi(1:end,sec))];
    bt = b(sec);
    Loca = 0:length(Vmode)-1;
    maxAbsValue = max(abs(Vmode));  % Find the maximum absolute value
    normalizedV(i,:) = Vmode;  % Divide each element by the maximum absolute value

    plot(1:r,abs(imag(Olam)),'.k')
    hold on
    end
end
%% Plot the mode shapes 
% surf(Loca,pt_T,normalizedV)
% caxis([0 0.7])
% shading interp
% view([0,90])



