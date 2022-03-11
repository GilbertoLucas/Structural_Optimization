using Optim, LinearAlgebra

# Augmented Lagrange Method code for inequality contrataint
function Augmented_Lagrange_Mult_Inequality_Box(f, g, lower, upper, x, k_max; α=1000., γ=1.01, tol=1e-1000)   
    λ, nite = zeros(length(g(x))), 0 
    for k in 1 : k_max; nite = k
        ψ(x) = max.(g(x),-λ./(2α))
        F(x) = f(x) + λ ⋅ g(x) + α * sum(ψ(x).^2)
        #x = Optim.minimizer(optimize(F, lower, upper, x, Fminbox(NelderMead()))) # Ordem zero
        x = Optim.minimizer(optimize(F, lower, upper, x, Fminbox(BFGS());autodiff = :forward)) # Ordem 1
        λ .+= 2*α*ψ(x)        
        if norm(ψ(x),Inf)<=tol; break; end # Ver Arora pg 481
        α *= γ;
    end
    return x, nite, α, λ
end;


#modulo de elasticidade
E = 3.0E7;
#altura da trelica
L = 10.0;
#angulo do carregamento e carregamento
#cenário 1
#theta = -pi/4;P = 40000.0;
#cenário 2
#theta = -pi/2;P = 30000.0;
#cenário 3
theta = -3*pi/4;P = 20000.0;

#limite inferior das variáveis de projeto
lb=[0.1, 0.1, 0.1, -10.0, -10.0, -10.0];
#limite superior das variáveis de projeto
ub=[6.1, 6, 6, 10.0, 10.0, 10.0];

#chute inicial
#x0 = [A1,A2,A3,x1,x2,x3]
x0 = [6.01,6.0,6.0,-5.0,0.0,5.0];

#Area das barras
A1(x) = x[1];
A2(x) = x[2];
A3(x) = x[3];

#Comprimento das barras
L1(x) = sqrt(L ^ 2 + x[4] ^ 2);
L2(x) = sqrt(L ^ 2 + x[5] ^ 2);
L3(x) = sqrt(L ^ 2 + x[6] ^ 2);


#Matriz de rigidez
k11(x) = E * (A1(x) * x[4] ^ 2 / L1(x) ^ 3 + A2(x) * x[5] ^ 2 / L2(x) ^ 3 + A3(x) * x[6] ^ 2 / L3(x) ^ 3);
k12(x) = E * (A1(x) * L * x[4] / L1(x) ^ 3 + A2(x) * L * x[5] / L2(x) ^ 3 + A3(x) * L * x[6] / L3(x) ^ 3);
k21(x) = k12(x);
k22(x) = E * (A1(x) * L ^ 2 / L1(x) ^ 3 + A2(x) * L ^ 2 / L2(x) ^ 3 + A3(x) * L ^ 2 / L3(x) ^ 3);

K(x) = [k11(x) k12(x); k21(x) k22(x)];

#Vetor de forcas
f1 = P*cos(theta);
f2 = P*sin(theta);
F = [f1,f2];

#Calculo dos deslocamentos
u(x) = K(x)\F;

#Calculo das tensoes
S1(x) = -E / L1(x) ^ 2 * (u(x)[2] * L + u(x)[1] * x[4]);
S2(x) = -E / L2(x) ^ 2 * (u(x)[2] * L + u(x)[1] * x[5]);
S3(x) = -E / L3(x) ^ 2 * (u(x)[2] * L + u(x)[1] * x[6]);

#Restricoes de desigualdades
g1(x) = abs(S1(x)) / 0.5000e4 - 0.1e1;
g2(x) = abs(S2(x)) / 0.20000e5 - 0.1e1;
g3(x) = abs(S3(x)) / 0.50000e4 - 0.1e1;
g(x) = [g1(x), g2(x), g3(x)];

#Funcao Objetivo - Volume das barras
V(x) = A1(x)*L1(x) + A2(x)*L2(x) + A3(x)*L3(x);

#Solve
@time xₑ, k, α, λ = Augmented_Lagrange_Mult_Inequality_Box(V, g, lb,ub, x0, 100)
println("x* : ", xₑ, ", iterations: ", k, ", α : ", α, ", λ : ", λ)