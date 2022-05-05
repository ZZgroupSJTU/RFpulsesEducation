function PulseFreqResponseAxes(pulse,freq_i,freq_f,N,M0,lin_phase,is_wrap, ax1, ax2, ax3, ax4);
% SYNTAX: [Mxy, Mz] = GetPulseFreqResponse(pulse,freq_i,freq_f,N,M0,lin_phase,is_wrap);
%
% The current function plots the frequency response of the input pulse from
% frequency freq_i to frequency freq_f, with N steps inbetween. It assumes
% all the spins start from the initial condition M0 which is supplied (and
% is assumed to be a 3x1 vector).
%
% Note that for such a plot, the gradient fields in the pulse object aren't
% used, only pulse.RFamp, pulse.RFphase and pulse.tp play a role.
%
% Optional parameters:
%     lin_phase - linear phase correction
%     is_wrap   - set to 1 for phase wrapping (interval [0,2*pi])

if nargin<6
    lin_phase = 0;
    is_wrap = 0;
end;

if nargin<7
    is_wrap = 0;
end;

% Step I: Create a spin structure
cs_vector = linspace(freq_i, freq_f, N);
for k=1:N
    spinstemp(k).r = [0; 0; 0];
    spinstemp(k).M = M0;
    spinstemp(k).cs = cs_vector(k);
end;

% Apply pulse
spinsout = ApplyPulse(spinstemp,pulse);

% Extract output magnetization
for k=1:N
    M(:,k) = spinsout(k).M;
end;

% Calculate Mxy, Mxy_phase, Mx, My (with lin_phase taken into account)
Mxy = sqrt(M(1,:).^2 + M(2,:).^2);
Mxy_phase = phase(M(1,:)+i*M(2,:)) + cs_vector*lin_phase;
if is_wrap==1
    Mxy_phase_wrapped = mod(Mxy_phase,2*pi)/pi;
else
    Mxy_phase_wrapped = Mxy_phase/pi;
end;
Mx = Mxy.*cos(Mxy_phase); 
My = Mxy.*sin(Mxy_phase); 


% Plot results
axes(ax1);
plot(cs_vector, M(3,:));
title('Mz');
xlabel('kHz');
ylabel('a.u.');

axes(ax2);
plot(cs_vector, Mxy);
title('|Mxy|');
xlabel('kHz');
ylabel('a.u.');

axes(ax3);
plot(cs_vector, Mxy_phase_wrapped);
title('Phase of Mxy');
xlabel('kHz');
ylabel('Units of \pi');

axes(ax4);
plot(cs_vector, Mx, 'r');
hold on
plot(cs_vector, My, 'b');
xlabel('kHz');
ylabel('a.u.');
legend('Mx','My');
title('Mx (red), My (blue)');
