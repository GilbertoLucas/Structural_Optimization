function [C,Ceq]=GetConstraints(x,L,P,E,theta)
%area das barras
A1 = x(1);
A2 = x(2);
A3 = x(3); 

%comprimento das barras
L1 = sqrt(L ^ 2 + x(4) ^ 2);
L2 = sqrt(L ^ 2 + x(5) ^ 2);
L3 = sqrt(L ^ 2 + x(6) ^ 2);

%calculo da matriz de rigidez
k11 = E * (A1 * x(4) ^ 2 / L1 ^ 3 + A2 * x(5) ^ 2 / L2 ^ 3 + A3 * x(6) ^ 2 / L3 ^ 3);
k12 = E * (A1 * L * x(4) / L1 ^ 3 + A2 * L * x(5) / L2 ^ 3 + A3 * L * x(6) / L3 ^ 3);
k21 = k12;
k22 = E * (A1 * L ^ 2 / L1 ^ 3 + A2 * L ^ 2 / L2 ^ 3 + A3 * L ^ 2 / L3 ^ 3);

K = [k11 k12; k21 k22];
    
%calculo do vetor de forca
f = [P*cos(theta) ;P*sin(theta)];

%calculo dos deslocamentos
u = linsolve(K,f);

%calculo das tensoes
S1 = -E / L1 ^ 2 * (u(2) * L + u(1) * x(4));
S2 = -E / L2 ^ 2 * (u(2) * L + u(1) * x(5));
S3 = -E / L3 ^ 2 * (u(2) * L + u(1) * x(6));

%restricoes de desigualdades
g1 = abs(S1) / 0.5000e4 - 0.1e1;
g2 = abs(S2) / 0.20000e5 - 0.1e1;
g3 = abs(S3) / 0.50000e4 - 0.1e1;
C=[g1,g2,g3];

%restricoes de igualdade
Ceq=[];