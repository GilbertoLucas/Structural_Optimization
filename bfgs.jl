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


function bfgs(f,grad_f,x0;tol=1e-6,maxitr=10000)

    #init data
    ite = 0;x=x0;conv=false;err=0;
    xh = zeros(length(x0),maxitr+1);xh[:,1] = x0;

    #init bfgs parameters
    m = length(x)
    D = Matrix(1.0I,m,m)
    g = grad_f(x)

    #iteractions
    while ite < maxitr
        #line search
        d = -D*g

        #Inexact line search (strongWolf)
        alpha = strong_backtracking(f,grad_f,x,d;α=1, β=1e-4, σ=0.1)

        #update x
        s = alpha*d;
        x += s
        
        gl = grad_f(x)
        y = gl - g
        D = D - (s*y'*D + D*y*s')/(s'*y) + (1+(y'*D*y)/(s'*y))[1]*(s*s')/(s'*y)

        #update bfgs parameters
        g = gl

        #update ite and x history
        ite += 1
        xh[:,ite+1] = x

        #local minimum point
        err = norm(s);
        if err < tol ; conv=true; break; end
    end
   
    println(ite+1)
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
#F1 = 0.0 ;F2 = 4.0; #nao ta convergindo

#Caso 3
F1 = 0.0; F2 = 400.0;

#function
F(u) = 0.5e0 * ka * (sqrt(u[1] ^ 2 + (La - u[2]) ^ 2) - La) ^ 2 + 0.5e0 * kb * (sqrt(u[1] ^ 2 + (Lb + u[2]) ^ 2) - Lb) ^ 2 - F1 * u[1] - F2 * u[2];
grad_F(u) = ForwardDiff.gradient(F,u)

#tolerancia
tol = 1e-8

#numero maximo de iteracoes
iterac_max = 10000

#ponto inicial
#u0 = [-0.01,0.0]
u0 = [ 0.00,0.0]

#u*
@time u_est = bfgs(F,grad_F,u0;tol=tol,maxitr=iterac_max)
u_est = u_est[1][:,end]
println("u*: ",u_est)
println("f(u*): ",F(u_est))

#gradiente no ponto de minimo
println("∇f(u*): ",grad_F(u_est))

#verificar a hessiana no ponto de minimo
hess_F(u) = ForwardDiff.hessian(F,u)
println("H(u*): ",hess_F(u_est))
println("λ: ",eigvals(hess_F(u_est)))