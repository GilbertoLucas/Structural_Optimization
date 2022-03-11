clc, clear, close all

%chute inicial
% x0 = [A1,A2,A3,x1,x2,x3]
x0 = [6,6,6,-5,0,5]';

%modulo de elasticidade
E = 3E7;
%altura da trelica
L = 10;

%angulo do carregamento e carregamento
%cenário 1
theta = -pi/4;P = 40000.0;
%cenário 2
%theta = -pi/2;P = 30000.0;
%cenário 3
%theta = -3*pi/4;P = 20000.0;

%restrições de desigualdades lineares
A=[];
b=[];
%restrições de igualdades lineares
Aeq=[];
beq=[];

%limite inferior das variáveis de projeto
lb=[0.1 0.1 0.1 -10.0 -10.0 -10.0]';
%limite superior das variáveis de projeto
ub=[6.1 6.1 6.1 10.0 10.0 10.0]';

%solve
options=optimset('Algorithm','active-set','Display','off');
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(@(x)GetVolume(x,L),x0,A,b,Aeq,beq,lb,ub,@(x)GetConstraints(x,L,P,E,theta),options);
[C,Ceq]=GetConstraints(x,L,P,E,theta);
if exitflag<=0, error('Problem did not converge - Check the ''exitflag'' fortroubleshooting'), end
fprintf('V=%.4f\n\tA1=%.4f\tA2=%.4f\tA3=%.4f\n\tx1=%.4f\tx2=%.4f\tx3=%.4f\n',fval,x);


%Grafico
%configuração  inicial
b1_x0 = [x0(4) 0]; b1_y = [L 0];
b2_x0 = [0 0];     b2_y = [L 0];
b3_x0 = [x0(6) 0]; b3_y = [L 0];

%configuração otimizada
b1_x = [x(4) 0];
b2_x = [x(5) 0];
b3_x = [x(6) 0];

plot(b1_x0,b1_y,'b',b2_x0,b2_y,'b',b3_x0,b3_y,'b',b1_x,b1_y,'r',b2_x,b2_y,'r',b3_x,b3_y,'r','linewidth',2);


