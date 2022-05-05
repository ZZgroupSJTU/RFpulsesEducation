function PlotOffsetTimeEvolution(pulse,offset,initState, position)
% Function name: PlotOffsetTimeEvolution
% SYNTAX: PlotOffsetTimeEvolution(pulse, offset, initState, position)
%
% Description: this function simulates the time evolution of a particular
% spin with a given offset (using the Bloch equations solver), and plots 
% the results: Mx(t), My(t), Mz(t), |Mxy(t)| and phase(Mxy(t)). 
% Relaxation is neglected.
%
% Variables Name     Type         Units       Description
% pulse              pulse        -           Input pulse structure
% offset             real 1x1     kHz         Offset of spin in question
% initState          real 3x1     -           Initial vector state [1]
% position           real 3x1     mm          Position of the spin [2]
%
% [1] State vectors are given as normalized column vectors. This is an
%     optional argument and is set to [0; 0; 1] if none is provided.
% [2] Position of the spin. Relevant if gradients are nonzero. Optional
%     argument; default value is [0; 0; 0].

if ~exist('initState','var')
    initState = [0; 0; 1];
end

if ~exist('position','var')
    position = [0; 0; 0];
end

% Apply pulse
[Mx, My, Mz] = ApplyPulseDiagnostics(offset, position, initState, pulse);

% Calculate timing axis
timeAxis = linspace(0, pulse.tp, length(pulse.RFamp));

% Plot results
figure

subplot(2,3,1);
plot(timeAxis, Mx);
title('Mx');
xlabel('ms');
ylabel('a.u.');

subplot(2,3,2);
plot(timeAxis, My);
title('My');
xlabel('ms');
ylabel('a.u.');

subplot(2,3,3);
plot(timeAxis, Mz);
title('Mz');
xlabel('ms');
ylabel('a.u.');

subplot(2,3,4);
plot(timeAxis, abs(Mx+i*My));
title('|M_{xy}|');
xlabel('ms');
ylabel('a.u.');

subplot(2,3,5);
plot(timeAxis, phase(Mx+i*My));
title('phase');
xlabel('ms');
ylabel('radians');