using ForwardDiff, Roots, Plots, LinearAlgebra, Optim

function strong_backtracking(f, ∇, x, d; α=1, β=0.1, σ=0.4)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN
    # bracket phase
    while true
        y = f(x + α*d)
        if y > y0 + β*α*g0 || (!isnan(y_prev) && y ≥ y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = ∇(x + α*d)⋅d
        if abs(g) ≤ -σ*g0
            return α
        elseif g ≥ 0
            αlo, αhi = α, α_prev
            break
        end
        y_prev, α_prev, α = y, α, 2α
    end
    
    # zoom phase
    ylo = f(x + αlo*d)
    while true
        α = (αlo + αhi)/2
        y = f(x + α*d)
        if y > y0 + β*α*g0 || y ≥ ylo
            αhi = α
        else
            g = ∇(x + α*d)⋅d
            if abs(g) ≤ -σ*g0
                return α
            elseif g*(αhi - αlo) ≥ 0
                αhi = αlo
            end
            αlo = α
        end
    end
end




f(x) = ((x[1]-1.0)^2) + (0.1*(x[2]-4)^2) + sin(5*pi*x[1])+sin(5*pi*x[2]) 
grad_f(x) = ForwardDiff.gradient(f,x)

# Ponto inicial
xk = [-1.0,1.0]

#direcao de busca
dk = -grad_f(xk)/norm(grad_f(xk))

beta = 0.1
alpha_f = 1.0

phi(alpha) = f(xk+alpha*dk)
d_phi(alpha) = grad_f(xk+alpha*dk)'dk

l_alpha(alpha) = phi(0)+alpha*beta*d_phi(0)
sigma = 0.4

wolf = []

for alpha = 0:0.01:alpha_f
    if (phi(alpha) <= l_alpha(alpha)) && (abs(grad_f(xk+alpha*dk)'dk) <= -sigma*grad_f(xk)'dk)
        push!(wolf,alpha)
    end
end


wolfbt = strong_backtracking(f,grad_f,xk,dk)

plot(phi,0,alpha_f,lw=2,lab="phi")
plot!(l_alpha,0,alpha_f,lw=2,lab="L_phi",c=:black)
scatter!(wolf,phi.(wolf),lab="StrongWolfe",ms=3,c=:red)
scatter!([wolfbt],[phi(wolfbt)],lab="StrongWolf - BT",c=:yellow)
plot!(legend=:outerbottomright)
plot!(title="beta = 0.1 e sigma = 0.4",xlabel = "alpha", ylabel = "phi(alpha), l_alpha(alpha")