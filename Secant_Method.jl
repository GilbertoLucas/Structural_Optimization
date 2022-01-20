using  LinearAlgebra
using Plots

function secant_method(f′, x0, x1, ϵ)
    g0 = f′(x0)
    Δ = Inf
    while abs(Δ) > ϵ
        g1 = f′(x1)
        Δ = (x1 - x0)/(g1 - g0)*g1
        x0, x1, g0 = x1, x1 - Δ, g1
    end
    return x1
end

f = x->sin(10x)+cos(3x)
df(x) = 10*cos(10*x)-3*sin(3*x)

#dois pontos iniciais
x0 = 3.4
x1 = 3.8

#convergencia
convergence = 1e-8      

xs = secant_method(df,x0,x1,convergence)
xs,f(xs),df(xs)

plot(f,0.,8.,lw=2,draw_arrow="true")
#scatter!([a,c],[f(a),f(c)])
scatter!([x0,x1],[f(x0),f(x1)])
scatter!([xs],[f(xs)])