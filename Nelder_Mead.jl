using LinearAlgebra, Statistics

function nelder_mead(f, S, ϵ; α=1.0, β=2.0, γ=0.5)
    Δ, y_arr = Inf, f.(S)
    iteractions = 0
    while Δ > ϵ
        p = sortperm(y_arr) # sort lowest to highest
        S, y_arr = S[p], y_arr[p]
        xl, yl = S[1], y_arr[1] # lowest
        xh, yh = S[end], y_arr[end] # highest
        xs, ys = S[end-1], y_arr[end-1] # second-highest
        xm = mean(S[1:end-1]) # centroid
        xr = xm + α*(xm - xh) # reflection point
        yr = f(xr)
    
        if yr < yl
            xe = xm + β*(xr-xm) # expansion point
            ye = f(xe)
            S[end],y_arr[end] = ye < yr ? (xe, ye) : (xr, yr)
        elseif yr > ys
            if yr ≤ yh
                xh, yh, S[end], y_arr[end] = xr, yr, xr, yr
            end
            xc = xm + γ*(xh - xm) # contraction point
            yc = f(xc)
            if yc > yh
                for i in 2 : length(y_arr)
                    S[i] = (S[i] + xl)/2
                    y_arr[i] = f(S[i])
                end
            else
                S[end], y_arr[end] = xc, yc
            end
        else
            S[end], y_arr[end] = xr, yr
        end

        Δ = std(y_arr, corrected=false)
        iteractions = iteractions + 1
    end
    println("Iteractions: ",iteractions)
    return S,iteractions
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

#tolerancia
tol = 1e-8

#ponto inicial
u0 = [-0.01,0.0]
#u0 = [0.0,0.0]
#numero de pontos iniciais = dimensao do problema + 1
u1 = [u0[1]-0.00001,u0[2]+0.00001]
u2 = [u0[1]+0.00001,u0[2]+0.00001]
u = [u0,u1,u2]

#function
f(u) = 0.5e0 * ka * (sqrt(u[1] ^ 2 + (La - u[2]) ^ 2) - La) ^ 2 + 0.5e0 * kb * (sqrt(u[1] ^ 2 + (Lb + u[2]) ^ 2) - Lb) ^ 2 - F1 * u[1] - F2 * u[2];

#tolerancia
tol = 1e-8

@time u_est,iterac = nelder_mead(f, u, tol)
println("u*: ",mean(u_est))
println("f(u*): ",f(mean(u_est)))
