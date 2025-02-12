function fhandle = plotCylinderX(VORT,s,i)

fhandle = figure (1)
subplot(1,s,i)

V2 = VORT;
A_rotated = rot90(V2);
imagesc(A_rotated)
