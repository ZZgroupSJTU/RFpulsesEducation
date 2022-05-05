function spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1)
% SYNTAX: 
%
%   spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1)
% 
% Creates a sample of J-coupled metabolites.
%         
% Input Variables
% Variable Name             Size          Units   Description
% csCenter                  1x1           ppm     Center frequency of Tx/Rx
% B0                        1x1           Tesla   Main field
% B1                        1x1                   B1-Scaling (global)
% isSecular                 1x1           -       Boolean. Set to 0 or 1 to signify
%                                                 if J coupling should be taken to be
%                                                 secular (only Ijz*Ikz terms) or not.
% linewidth                 1x1           Hz      An ad-hoc term which "relaxes" the
%                                                 acquired FID and sets the linewidth.
%                                                 (Functions as an "effective T2")
%
% Output Variables
% spins                                           Spin structure (with no molecules)
%                                     
% The spins structure is in the following format:
%
%   spins.csCenter                        ppm     Center of irradiation & acquisition
%   spins.B0                              Tesla   Main field strength
%   spins.B1                              [scale] B1 scaling factor (global)
%   spins.isSecular                       0,1     If 1, the J-coupling Hamiltonian 
%                                                 will only contain Ikz*Ijz terms.
%   spins.linewidth                       Hz      An effective T2 line broadening term
%   spins.molecule(i).csVec               ppm     The molecule's spins' chemical shifts
%   spins.molecule(i).JMatrix             Hz      The molecule's spins' J coupling
%                                                 matrix
%   spins.molecule(i).nucleus             -       Cell array containing
%                                                 spins' nuclei.
%   spins.molecule(i).spin(j).rho         -       The molecule's density matrix
%   spins.molecule(i).spin(j).r           mm      The spin's position (1x3 vector).
%   spins.molecule(i).spin(j).B1          [scale] B1 scaling factor (per molecule).
%   spins.molecule(i).spin(j).B0          ppm     Spin's local offset.
%   spins.molecule(i).spin(j).linewidth   Hz      Spin's own linewidth.

if nargin<1, csCenter = 0; end
if nargin<2, B0 = 3; end
if nargin<3, B1 = 1; end
if nargin<4, isSecular = 0; end
if nargin<5, linewidth = 5; end

spins.csCenter = csCenter;
spins.B0 = B0;
spins.B1 = B1;
spins.isSecular = isSecular;
spins.linewidth = linewidth;
spins.molecule = [];

