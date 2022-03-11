function [C,Ceq]=GetConstraints(x,B,P,Sy,r1,E)
r2=x(1);
H =x(2);


g2=P*sqrt(B^2+H^2)/(H*pi*r2^2*Sy)-1;
g3=4*P*(B^2+H^2)^(3/2)/(H*pi^3*r2^4*E)-1;

C=[g2 g3];


Ceq=[];

P*B/(pi*H*(r1^2)*Sy)