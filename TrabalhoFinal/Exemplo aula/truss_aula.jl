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

Sy=350.0;
r1=0.4;
E=2000.0;
P=200.0;
B=1.0;

f(x) = pi*(x[1]^2)*(x[2]^2+B^2)^(1/2);

g1(x) = P*B/(pi*x[2]*(r1^2)*Sy) - 1;
g2(x) = P*((B^2+x[2]^2)^(1/2))/(x[2]*pi*x[1]^2*Sy)-1;
g3(x) = 4*P*(B^2+x[2]^2)^(3/2)/(x[2]*pi^3*x[1]^4*E)-1;

g(x) = [g1(x), g2(x), g3(x)];

lb = [0.4, 1.0];
ub = [0.8, 6.0];

x0 = [0.5, 1.1]

@time xₑ, k, α, λ = Augmented_Lagrange_Mult_Inequality_Box(f, g, lb,ub, x0, 100)
println("x* : ", xₑ, ", iterations: ", k, ", α : ", α, ", λ : ", λ)