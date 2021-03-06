function SLRPlotMagnetizationFromAB(A,B, initMag, numPoints)
% Given the two SLR polynomials, A and B, this computes the final
% magnetization (given initial magnetization) as a function of the
% normalized frequency, and then plots its Mz, Mxy components.


freqAxis = [-1/2:1/numPoints:1/2 - 1/numPoints];
z = exp(1i*2*pi*freqAxis);

AN = 0.*freqAxis;
BN = 0.*freqAxis;

for k=1:length(A)
    AN = AN + A(k)*z.^(k-1);
    BN = BN + B(k)*z.^(k-1);
end

a = AN.*z.^(length(A)/2);
b = BN.*z.^(length(A)/2);


Q{numPoints}=zeros(2,2);
for k=1:length(freqAxis)
    Q{k} = [a(k) -conj(b(k)); 
            b(k)  conj(a(k))];
    P = Q{k}*[1; 0];

    PPP=Q{k}*(Q{k}');
    fprintf('Step %d, PPP = %f \n', k, PPP(1,1));

    % Convert spinor to 3D real magnetization
    x = P(1);
    y = P(2);
    M = [x*conj(y) + conj(x)*y;
            1i*(conj(x)*y - x*conj(y));
            abs(x)^2 - abs(y)^2];
    Mxy(k) = M(1)+1i*M(2);
    Mx(k) = M(1);
    My(k) = M(2);
    Mz(k) = M(3);
end

figure

subplot(2,2,1)
plot(freqAxis, Mx);
hold
plot(freqAxis, My,'r');
title('Mx, My');

subplot(2,2,2)
plot(freqAxis, My);
title('My');

subplot(2,2,3)
plot(freqAxis, abs(Mxy));
title('Mxy');

subplot(2,2,4)
plot(freqAxis, phase(Mxy));
title('Phase(Mxy)');

