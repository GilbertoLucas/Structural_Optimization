using Plots

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
            return a < c ? (a,b, c) : (c, b, a)
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


#definicao da funcao
f = x-> -(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x    #multiplicada por -1 para encontrar o minimo
g = x-> -(-(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x) #funcao original

#definiciao do ponto inicial
x0 = 0
#busca do intervalo
a,b,c = bracket_minimum(f,x0)

#numero de iteracoes
φ = MathConstants.φ
n=(c-a)/(10^-3*log(φ))

am,bm,cm = quadratic_fit_search(f,a,b,c,n)
#(am,f(am))
#am,bm,cm

println("x* = ", cm)

plot(g,a-1.,c+1.,lw=2,draw_arrow="true",label="Function")
scatter!([a,b,c],[g(a),g(b),g(c)],label="interval")
scatter!([cm],[g(cm)],label="x*")
