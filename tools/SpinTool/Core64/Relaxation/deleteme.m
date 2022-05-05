
spins = InitSpins3D('numSpins', [1 1 1], 'initMag', [0 0 1], 'B1', 1.0);
seq = {{'hard', 180, 0, 1, 'normal', 1}};
spinsOut = ApplySequence3D(spins, seq);
squeeze(spinsOut.M)
