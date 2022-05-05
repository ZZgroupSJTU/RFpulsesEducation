function foo = b1rf2varian(rf_fm,rf_am,dt,fname,nramp)

% the format goes as: phase (deg) \t amp (0-1023) (must be positive? 
% : write 2 files, one with all positive + a pi phase shift, and one where it can go negative)
% \t relative width. one decimal place for each

% integrate the fm waveform to get phase
phs_base = cumsum(rf_fm(:)*dt*360); % degrees

ramp = [0:(nramp+2)/(nramp+2-1):nramp+2]'./(nramp+2);
ramp = ramp(2:end-1);

rf_amr = [-ramp;ones(nramp/2,1);rf_am(:);ones(nramp/2,1);-flipud(ramp)];

rf_phsr = [zeros(nramp*1.5,1);phs_base(:);zeros(nramp*1.5,1)];
rf_phsr = rf_phsr + 180*(rf_amr < 0);
rf_phsr(rf_phsr < 0) = rf_phsr(rf_phsr < 0)+360;

foo = rf_amr; % convenience for plotting the am waveform

rf_amr = abs(rf_amr);

disp(sprintf('Pulse duration (us): %f',round(dt*1000000*length(rf_amr))));

% open file
fp = fopen([fname '.RF'],'w');

% write header
fprintf(fp,'# @(#)%s.RF\n',fname);
fprintf(fp,'# \n');
fprintf(fp,'# TB: A B1-selective pulse\n');
fprintf(fp,'# ***************************************************\n');
fprintf(fp,'# %s.RF\n',fname);
fprintf(fp,'# VERSION\t1.0\n');
fprintf(fp,'# TYPE\t\tselective\n');
fprintf(fp,'# MODULATION\tamplitude\n');
fprintf(fp,'# EXCITEWIDTH\t1.0\n');
fprintf(fp,'# INVERTWIDTH\t1.0\n');
fprintf(fp,'# INTEGRAL\t%0.4f\n',sum(rf_amr)/length(rf_amr));
fprintf(fp,'# ***************************************************\n');

for ii = 1:length(rf_amr)
  fprintf(fp,'%0.1f\t%0.1f\t%0.1f\n',rf_phsr(ii),rf_amr(ii)/max(abs(rf_amr))*1023,1);
end

fclose(fp);

%eval(sprintf('!scp %s.RF grissowa@nmr47:~/vnmrsys/shapelib/',fname));

