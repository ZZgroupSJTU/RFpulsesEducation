% Single spin operators
Ix = [0 1; 1 0]/2;
Iy = [0 -i; i 0]/2;
Iz = [1 0; 0 -1]/2;
I = eye(2);


% Note about Hamiltonians:
%
% Single Spin Hamiltonian
% -----------------------
% H = M dot B
% For RF = 0, Bz = B0, we have:
% H = Mz*Bz = gamma*Sz*B0 = gamma*hbar*sigma_z*B0  
%
% J-Coupling Hamiltonian
% ----------------------
% H = M1 dot B1 + M2 dot B2 + 2*pi*J*(S1 dot S2)
