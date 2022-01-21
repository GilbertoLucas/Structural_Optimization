using ForwardDiff, Roots, Plots, LinearAlgebra, Optim

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


armj = []

for alpha = 0:0.01:alpha_f
    if phi(alpha) <= l_alpha(alpha)
        push!(armj,alpha)
    end
end

plot(phi,0,alpha_f,lw=2,lab="phi")
plot!(l_alpha,0,alpha_f,lw=2,lab="L_phi",c=:black)
scatter!(armj,phi.(armj),lab="Armijo",ms=3,c=:red)
plot!(legend=:outerbottomright)
plot!(title="beta = 0.1",xlabel = "alpha", ylabel = "phi(alpha), l_alpha(alpha")