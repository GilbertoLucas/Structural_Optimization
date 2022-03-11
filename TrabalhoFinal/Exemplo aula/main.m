clc, clear, close all

Sy=350;
r1=0.5;
E=2000;
P=200;
B=1;


A=[0 -pi*r1^2*Sy/(P*B)];
b=[-1];
Aeq=[];
beq=[];
lb=[0.4 1]';
ub=[0.8 6]';

x0=[0.5 1]';

options=optimset('Algorithm','active-set','Display','off');
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(@(x)GetVolume(x,B),x0,A,b,Aeq,beq,lb,ub,@(x)GetConstraints(x,B,P,Sy,r1,E),options);
[C,Ceq]=GetConstraints(x,B,P,Sy,r1,E);
if exitflag<=0, error('Problem did not converge - Check the ''exitflag'' fortroubleshooting'), end
fprintf('V=%.4f\tr2=%.4f\tH=%.4f\n',fval,x);
fprintf('Ax-b=%.4f\tnonlcon=[%.4f %.4f]\n',A*x-b,C);