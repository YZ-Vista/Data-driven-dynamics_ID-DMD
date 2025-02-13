function fhandle = plotCylinderX_Denoise(VORT,s,i)

fhandle = figure (1)
subplot(1,s,i)
V2 = VORT;

vortmin = 0.1;
% % normalize values... not symmetric
V2(V2<vortmin) = vortmin;
A_rotated = flip(V2,1);
imagesc(A_rotated)
