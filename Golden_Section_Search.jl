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
    


f = x->sin(10x)+cos(3x)
x0 = 3.4
#busca do intervalo
a,c = bracket_minimum(f,x0)

#numero de iteracoes
n=(c-a)/(10^-3*log(φ))

am,bm = golden_section_search(f,a,c,n)
(am,f(am))
am,bm


plot(f,0.,4.,lw=2,draw_arrow="true")
scatter!([a,c],[f(a),f(c)])
scatter!([am,bm],[f(am),f(bm)])
