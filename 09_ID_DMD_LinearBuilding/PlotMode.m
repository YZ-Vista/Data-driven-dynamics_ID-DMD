function PlotMode(BUILD)

Vmode = [0; BUILD];
Loca = 0:length(Vmode)-1;

maxAbsValue = max(abs(Vmode));  % Find the maximum absolute value
normalizedV = Vmode / maxAbsValue;  % Divide each element by the maximum absolute value

plot(normalizedV,Loca,'o--','LineWidth',1.2)
hold on
plot(0*normalizedV,Loca,'o-k','LineWidth',1.2)
hold on

xlabel('mode')
ylabel('floor')

yticks(Loca)
xlim([-1,1])