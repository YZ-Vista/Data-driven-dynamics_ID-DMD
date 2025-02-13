function fhandle = plotCylinderX(VORT,s,i)

fhandle = figure (1)
subplot(1,s,i)
V2 = VORT;

vortmin = 0;
vortmax = 100;

V2(V2<vortmin) = vortmin;
V2(V2>vortmax) = vortmax;
imagesc(V2)

