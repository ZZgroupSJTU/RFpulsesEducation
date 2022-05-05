% figure
% PlotBloch(0, 0,  [4.6217   18.2217    3.1003], 0);
% % colormap('jet');
% hold on
% xlim([-1 1]);
% B = [0.3 0 1.1];
% arrowSize = 0.017/norm(B);
% % B = B./norm(B);
% dt = 0.01;
% N = 40;
% T = 2*pi;
% dt = T/N;
% Mi = [0; 0; 1];
% C = [0.1 0.1 0.5];
% scaleFactor = 0.8;
% for idx=1:N
%     Mf = RotMat(B, 2*pi/N)*Mi;
%     X = [0 Mi(1) Mf(1)]*scaleFactor;
%     Y = [0 Mi(2) Mf(2)]*scaleFactor;
%     Z = [0 Mi(3) Mf(3)]*scaleFactor;
%     patch(X, Y, Z, C, 'EdgeColor', C, 'FaceAlpha', 0.5);
%     % plot3(X, Y, Z);
%     Mi = Mf;
% end
% h=arrowPlot([0 0 0], B*1.0+[0 0.01 0.01], [1 1 0], arrowSize);
% 
% set(h, 'FaceColor', [1 0 0]);


figure
PlotBloch(0, 0,  [ 13.1556    9.9475    7.3914], 0);
% colormap('jet');
hold on
xlim([-1 1]);
B = [0.3 0 1.1];
arrowSize = 0.017/norm(B);
% B = B./norm(B);
dt = 0.13;
N = 40;
% T = 1;
da = norm(B)*dt;
Mi = [0; 0; 1];
C = [0.1 0.1 0.5];
scaleFactor = 0.8;
for idx=1:N
    Mf = RotMat(B, da)*Mi;
    X = [0 Mi(1) Mf(1)]*scaleFactor;
    Y = [0 Mi(2) Mf(2)]*scaleFactor;
    Z = [0 Mi(3) Mf(3)]*scaleFactor;
    patch(X, Y, Z, C, 'EdgeColor', C, 'FaceAlpha', 0.5);
    % plot3(X, Y, Z);
    Mi = Mf;
end
h=arrowPlot([0 0 0], B*1.0+[0 0.01 0.01], [1 1 0], arrowSize);

set(h, 'FaceColor', [1 0 0]);


