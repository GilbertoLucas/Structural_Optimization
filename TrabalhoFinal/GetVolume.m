function V=GetVolume(x,L)
%area das barras
A1 = x(1);
A2 = x(2);
A3 = x(3); 

%comprimento das barras
L1 = sqrt(L ^ 2 + x(4) ^ 2);
L2 = sqrt(L ^ 2 + x(5) ^ 2);
L3 = sqrt(L ^ 2 + x(6) ^ 2);

%volume das barras (função objetivo)
V = A1*L1 + A2*L2 + A3*L3;