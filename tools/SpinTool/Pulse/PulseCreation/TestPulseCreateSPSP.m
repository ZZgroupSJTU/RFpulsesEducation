
Nc = 16; % Number of cycles
Tc = 2; % ms
Np = 64; % Number of points per cycle
G = 0.5; % kHz/mm

GCycle = [ones(1,Np/2)*G, -ones(1,Np/2)*G];

GShape = [-ones(1,Np/4)*G, repmat(GCycle, 1, Nc)];
N = numel(GShape);

flipAngle = 90;
numLobes = 4;
totTime = Tc*Nc; 
pulseFreq = PulseCreateSinc(numLobes, totTime, Nc, flipAngle);

flipAngle = 3;
numLobes = 4;
pulseSpat = PulseCreateSinc(numLobes, Tc/2, Np/2, flipAngle);


RFamp = zeros(1, Np/4);
RFphase = zeros(1, Np/4);

for idx=1:Nc
    RFamp = [RFamp, pulseSpat.RFamp.*pulseFreq.RFamp(idx), zeros(1,Np/2)];
    RFphase = [RFphase, pulseSpat.RFphase + pulseFreq.RFphase(idx), zeros(1,Np/2)];
end


p.tp = totTime; 
p.RFamp = RFamp*200;
p.RFphase = RFphase;
p.Gx = zeros(1, N);
p.Gy = zeros(1, N);
p.Gz = GShape;

freqLimits = [-0.3 0.3 61];
spatialLimits = [-20 20 71];
initMag = [0;1;0];
csVec = linspace(freqLimits(1), freqLimits(2), freqLimits(3));
posVec = linspace(spatialLimits(1), spatialLimits(2), spatialLimits(3));
purgeMoment = -45

counter = 0;
for idxFreq=1:freqLimits(3)
    for idxPos=1:spatialLimits(3)
        counter = counter + 1;
        spins(counter).r = [0; 0; posVec(idxPos)];
        spins(counter).M  = initMag;
        spins(counter).cs = csVec(idxFreq);  % in kHz!
        spins(counter).T1 = 1e6; % ms
        spins(counter).T2 = 1e6; % ms
        spins(counter).M0 = 1; % a.u.
        spins(counter).B1 = 1; % Scales RF
        spins(counter).B0 = 0; % Offset, in kHz
        spins(counter).RS = 1; % Receiver sensitivity
    end
end

spinsOut = ApplyPulseRelax(spins, p);
spinsOut = PurgeMoment(spinsOut, 0, 0, purgeMoment, 0.001);

counter = 0;
for idxFreq=1:freqLimits(3)
    for idxPos=1:spatialLimits(3)
        counter = counter + 1;
        Mx(idxFreq, idxPos)  = spinsOut(counter).M(1);
        My(idxFreq, idxPos)  = spinsOut(counter).M(2);
        Mz(idxFreq, idxPos)  = spinsOut(counter).M(3);
        Mxy(idxFreq, idxPos) = spinsOut(counter).M(1) + 1i*spinsOut(counter).M(2);
    end
end


figure
imagesc(posVec, csVec, Mx)
xlabel('mm');
ylabel('kHz');

%PlotPulseSpatialSpectralResponse(p, [0;0;1], [-0.30 0.30 61], [-20 20 71], {'mxy','mx'}, 'z', 'mesh');
