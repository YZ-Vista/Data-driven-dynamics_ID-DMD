function PlotMode(Mode)

Vmode = [0; Mode];
Loca = 0:length(Vmode)-1;

maxAbsValue = max(abs(Vmode));  % Find the maximum absolute value
normalizedV = Vmode/ maxAbsValue;  % Divide each element by the maximum absolute value

plot(Loca,normalizedV,'o--','LineWidth',1.2)
hold on
plot(Loca,0*normalizedV,'o-k','LineWidth',1.2)
hold on

xlabel('freq')
ylabel('mode')

xticks(Loca)
% ylim([-1,1])