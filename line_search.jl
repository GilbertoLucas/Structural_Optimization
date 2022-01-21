
function bracket_minimum(f, x=0; s=1e-2, k=2.0)
    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end
    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a,b, c) : (c,b, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
end

function quadratic_fit_search(f, a, b, c, n)
    ya, yb, yc = f(a), f(b), f(c)
    for i in 1:n-3
        x = 0.5*(ya*(b^2-c^2)+yb*(c^2-a^2)+yc*(a^2-b^2)) /
            (ya*(b-c) +yb*(c-a) +yc*(a-b))
        yx = f(x)
        if x > b
            if yx > yb
                c, yc = x, yx
            else
                a, ya, b, yb = b, yb, x, yx
            end
        elseif x < b
            if yx > yb
                a, ya = x, yx
            else
                c, yc, b, yb = b, yb, x, yx
            end
        end
    end
    return (a, b, c)
end

φ = MathConstants.φ
function golden_section_search(f, a, b, n)
    ρ = φ-1
    d = ρ * b + (1 - ρ)*a
    yd = f(d)
    for i = 1 : n-1
        c = ρ*a + (1 - ρ)*b
        yc = f(c)
        if yc < yd
            b, d, yd = d, c, yc
        else
            a, b = b, c
        end
    end
    return a < b ? (a, b) : (b, a)
end
    
function line_search(f, x, d)
    objective = α -> f(x + α*d)
    a, b,c = bracket_minimum(objective)
    
    #numero de iteracoes
    n=(c-a)/(10^-3*log(φ))
    
    #pelo metodo golden section search
    α,beta = golden_section_search(objective, a, c,n)

    #pelo metodo quadratic fit search
    #α,beta,gamma = quadratic_fit_search(objective, a, b, c,n)
    return (x + α*d,α)
end

#funcao multivariavel
f(x) = sin(x[1]*x[2])+exp(x[2]+x[3])-x[3]

#ponto inicial
x0 = [1,2,3]

#direcao de busca
d = [0,-1,-1]

(x,α) = line_search(f, x0, d)
println("x* = ", x)