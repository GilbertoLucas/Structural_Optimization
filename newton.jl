using ForwardDiff, Roots, Plots, LinearAlgebra, Optim,Statistics

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

function newton(f,grad_f,H,x0,tol=1e-6,maxitr=100)
    
    #init data
    ite = 0;x=x0;conv=false;err=0;
    xh = zeros(length(x0),maxitr+1);xh[:,1] = x0

    #iteractions
    while ite < maxitr
        
        #search direction
        g = grad_f(x)
        d = -H(x)\g
        if any(isfinite.(d).==false);break;end

        #line search
        alpha = strong_backtracking(f,grad_f,x,d;α=1, β=1e-4, σ=0.4)

        #update x 
        delta_x = alpha*d
        x += delta_x

        #update ite and x history
        ite += 1
        xh[:,ite+1] = x

        #local minimum point
        err = norm(delta_x)
        if err < tol ; conv=true; break; end
    end

    println("Iteractions: ", ite)
    #return x history. convergence and error
    return xh[:,1:ite+1], conv, err
end

#valores dos parametros
ka = 9.0 
kb = 2.0
La = 10.0
Lb = 10.0

#Caso 1
#F1 = 2.0; F2 = 4.0;

#Caso 2
#F1 = 0.0 ;F2 = 4.0;

#Caso 3
F1 = 0.0; F2 = 400.0;

#function
F(u) = 0.5e0 * ka * (sqrt(u[1] ^ 2 + (La - u[2]) ^ 2) - La) ^ 2 + 0.5e0 * kb * (sqrt(u[1] ^ 2 + (Lb + u[2]) ^ 2) - Lb) ^ 2 - F1 * u[1] - F2 * u[2];
grad_F(u) = ForwardDiff.gradient(F,u)
Hess_F(u) = ForwardDiff.hessian(F,u)

#tolerancia
tol = 1e-8

#numero maximo de iteracoes
iterac_max = 200

#ponto inicial
u0 = [-0.01,0.0]
#u0 = [ 0.00,0.0]

@time u_est = newton(F,grad_F,Hess_F,u0,tol,iterac_max)
u_est= u_est[1][:,end]
println("u*: ",u_est)
println("f(u*): ",F(u_est))

#gradiente no ponto de minimo
println("∇f(u*): ",grad_F(u_est))

#verificar a hessiana no ponto de minimo
println("H(u*): ",Hess_F(u_est))
println("λ: ",eigvals(Hess_F(u_est)))