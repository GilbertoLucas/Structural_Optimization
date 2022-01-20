

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

f = x->sin(10x)+cos(3x)
x0 = 3.4
#busca do intervalo
a,c = bracket_minimum(f,x0)

using Plots

plot(f,0.,4.,lw=2,draw_arrow="true")
scatter!([a,c],[f(a),f(c)])