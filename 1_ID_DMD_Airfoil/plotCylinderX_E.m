function fhandle = plotCylinderX_E(VORT,i)

fhandle = figure (1)
subplot(1,3,i)
V2 = VORT;

vortmin = 0;
vortmax = 0.1;

V2(V2<vortmin) = vortmin;
V2(V2>vortmax) = vortmax;
A_rotated = flip(V2,1);
imagesc(A_rotated)

