function [theta, phi] = inverseSLR(a,b)
% InverseSLR uses the inverse SLR algorithm to determine the series of 
% hard pulses (flip angle and phase) that produce the polynomials A and B
% Adapted from the thesis of P. Le Roux (algorithm 2, p. 62)


% Initialise the arrays that hold the c and s functions
n = length(b);
c = zeros(1,n);
s = zeros(1,n);

% Iterate from the last pulse to the first
for k = n:-1:1

    % The value of c and s can be obtained from the coefficients 
    % of A and B of zeroth order  
    norm = sqrt( conj(b(1))*b(1) + conj(a(1))*a(1) );
    cc = a(1)/norm;
    ss = b(1)/norm;
    
    % The recursion relations provide A(i-1) and B(i-1) 
    aprev = conj(cc) * a + conj(ss) * b;
    b = -ss*a+cc*b;
    a(1:k-1) = aprev(1:k-1);
    b(1:(k-1)) = b(2:k);
    
    % The 
    c(k) = cc;
    s(k) = ss;
    
end

% c and s can be converted into theta and phi
theta = 2*atan2(abs(s),abs(c));
phi = angle(s./c);

end