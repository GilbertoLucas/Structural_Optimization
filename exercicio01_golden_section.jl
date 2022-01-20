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
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
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
    
#definicao da funcao
f = x-> -(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x    #multiplicada por -1 para encontrar o minimo
g = x-> -(-(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x) #funcao original

#ponto inicial
initial_point = 0
#busca do intervalo
a,c = bracket_minimum(f,initial_point)

#numero de iteracoes
n=(c-a)/(10^-3*log(φ))

am,bm = golden_section_search(f,a,c,n)
#(am,f(am))
#am,bm

println("x* = ", am)

plot(g,a-1.,c+1.,lw=2,draw_arrow="true",label="Function")
scatter!([a,c],[g(a),g(c)],label="interval")
scatter!([am,bm],[g(am),g(bm)],label="x*")


