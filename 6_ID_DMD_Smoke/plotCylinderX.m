function fhandle = plotCylinderX(VORT,i)

fhandle = figure (1)
subplot(1,3,i)
V2 = VORT;

A_rotated = flip(V2,1);
imagesc(A_rotated)
