%% generate figure for lipid quadratic model

CNgg = [0, 12,24,30,36,48,60,90,120,180,240];
lipid = [0, 13.71, 13.48, 19.15, 15.98, 24.09, 33.78, 38.09 44.36, 41.66, 26.64];
error_bar = [0, 1.10, 1.90, 1.53, 3.4, 0.8, 4.68, 8.28, 2.56, 5.74, 3.21];
f = shadedErrorBar(CNgg,lipid,error_bar,'lineProps','-g');

f.mainLine.Color = [77/255 190/255 238/255];
f.mainLine.LineStyle = ':' ;
f.mainLine.Marker = 'x' ;
f.mainLine.MarkerSize = 20;
f.mainLine.LineWidth = 4;
f.edge = 'No'; % adjust manually
f.patch.FaceColor = [77/255 190/255 238/255];
% add equation from editPlot 
legend boxoff


ylabel('% Lipid content in DCW (w/w)')
xlabel('C/N (g/g)')
set(gca,'FontSize',50)
xlim([0 240])