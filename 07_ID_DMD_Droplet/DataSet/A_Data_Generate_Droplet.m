clear all, 
clc
load a15V_1050_150_100.mat imageData;
Droplet{1} = imageData; Vol{1} = 15;
load a18V_1050_150_100.mat imageData;
Droplet{2} = imageData; Vol{2} = 18;
load a21V_1050_150_100.mat imageData;
Droplet{3} = imageData; Vol{3} = 21;
load a24V_1050_150_100.mat imageData;
Droplet{4} = imageData; Vol{4} = 24;
load a27V_1050_150_100.mat imageData;
Droplet{5} = imageData; Vol{5} = 27;
load a30V_1050_150_100.mat imageData;
Droplet{6} = imageData; Vol{6} = 30;

imagesc(Droplet{1}(:,:,1))
[lx,ly,td] = size(Droplet{1});
for s = 1:6
    for i = 1:td
         Xdroplet{s}(:,i) = reshape(Droplet{s}(:,:,i),lx*ly,1);
    end
end
save DropLet_D.mat Xdroplet Vol lx ly
