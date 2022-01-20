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

f = x->sin(10x)+cos(3x)
df(x) = 10*cos(10*x)-3*sin(3*x)
ddf(x) = -100*sin(10x)-9*cos(3x)

#converge para o ponto de minimo
initial_point = 5.1       #ponto inicial
convergence = 1e-8      #convergencia
maximum_iterations = 10 #numero maximo de iteracoes
xs = newtons_method(df,ddf,initial_point,convergence,maximum_iterations)
xs,f(xs),ddf(xs),df(xs)


#converge para o ponto de maximo
#initial_point = 4.45       #ponto inicial
#convergence = 1e-8      #convergencia
#maximum_iterations = 10 #numero maximo de iteracoes
#xs = newtons_method(df,ddf,initial_point,convergence,maximum_iterations)
#xs,f(xs),ddf(xs),df(xs)

plot(f,0.,8.,lw=2,draw_arrow="true")
#scatter!([a,c],[f(a),f(c)])
scatter!([initial_point],[f(initial_point)])
scatter!([xs],[f(xs)])