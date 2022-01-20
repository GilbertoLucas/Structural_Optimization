using  LinearAlgebra
using Plots


function newtons_method(∇f, H, x, ϵ, k_max)
    k, Δ = 1, Inf
    while abs(Δ) > ϵ && k ≤ k_max
        Δ =  ∇f(x)/abs(H(x))
        x -= Δ
        k += 1
    end
    return x
end


#definicao da funcao
f = x-> -(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x    #multiplicada por -1 para encontrar o minimo
g = x-> -(-(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x) #funcao original

#definicao da primeira derivada
df(x) = 1.2*x/((x^2+1)^(3/2))+(1-0.4/(x^2+1))*x/((x^2+1)^(1/2))-1

#definicao da segunda derivada
ddf(x) = (-0.6*x^2+1.8)/((x^2+1)^(5/2))

#converge para o ponto de minimo
initial_point = 1       #ponto inicial
convergence = 1e-8      #convergencia
maximum_iterations = 10 #numero maximo de iteracoes
xs = newtons_method(df,ddf,initial_point,convergence,maximum_iterations)
#xs,f(xs),ddf(xs),df(xs)

println("initial_point = ", initial_point)
println("x* = ", xs)
println("f(x*) = ", g(xs))


plot(g,min(initial_point,xs)-1.,max(initial_point,xs)+1.,lw=2,draw_arrow="true",label="Function")
scatter!([initial_point],[g(initial_point)])
scatter!([xs],[g(xs)])
