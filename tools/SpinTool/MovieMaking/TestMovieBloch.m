N=200;
%pulse = PulseCreateSinc(4, 8, 128, 90);
pulse = PulseCreateConst(1, N, 3, 0);
%pulse = [zeros(1,N); zeros(1,N); ones(1,N)]*(5);
pulseParameter = 0.9;
spinPos = 0;
initialMag = [1/sqrt(2); 0; 1/sqrt(2)];
numFramesToSkip = 1;
addAnnotations = 0;
plotJustRF = 0;
isPlotPlaneRF = 0;
cameraPosition = [0.2; 0.5; 0.5];

outputMovie=MovieBloch(pulse,pulseParameter,spinPos,initialMag,numFramesToSkip,addAnnotations, plotJustRF, isPlotPlaneRF);

movie(outputMovie);

return

vidObj = VideoWriter('C:\Users\Assaf\Desktop\Gradient5.avi');
open(vidObj);
for k=1:200, writeVideo(vidObj, outputMovie(k)); end; 
close(vidObj)
