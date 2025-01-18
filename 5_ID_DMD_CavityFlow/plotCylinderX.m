function fhandle = plotCylinderX(VORT,i)

fhandle = figure (1)
subplot(1,3,i)

V2 = VORT;

A_rotated = rot90(V2);
imagesc(A_rotated)
