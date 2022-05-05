function spinsOut = ConcatSpins(spins1, spins2)

spinsOut = spins1;

spinsOut(numel(spins1)+numel(spins2)).r=[0; 0; 0];
for idx=1:numel(spins2)
    spinsOut(numel(spins1)+idx)=spins2(idx);
end
