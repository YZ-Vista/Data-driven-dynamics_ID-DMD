load Bur_Train.mat x tT Burg_train Param
parameter = Param{1} % Burg P*t: row: Position; Column: Time

surf(tT,x,Burg_train{1}), 
view([90,-90])
colorbar