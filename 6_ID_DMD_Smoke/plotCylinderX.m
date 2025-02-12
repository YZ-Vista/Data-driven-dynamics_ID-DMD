function fhandle = plotCylinderX(VORT,s,i)

fhandle = figure (1)
subplot(1,s,i)
V2 = VORT;

A_rotated = flip(V2,1);
imagesc(A_rotated)
