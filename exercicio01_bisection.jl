using Plots

function bracket_sign_change(f′, a, b; k=2)
    if a > b; a,b = b,a; end # ensure a < b
    center, half_width = (b+a)/2, (b-a)/2
    while f′(a)*f′(b) > 0
        half_width *= k
        a = center - half_width
        b = center + half_width
    end
    return (a,b)
end

function bisection(f′, a, b, ϵ)
    if a > b; a,b = b,a; end # ensure a < b
    
    ya, yb = f′(a), f′(b)
    if ya == 0; b = a; end
    if yb == 0; a = b; end
    
    while b - a > ϵ
        x = (a+b)/2
        y = f′(x)
        if y == 0
            a, b = x, x
        elseif sign(y) == sign(ya)
            a = x
        else
            b = x
        end
    end
    return (a,b)
end
    

    
#definicao da funcao
f = x-> -(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x    #multiplicada por -1 para encontrar o minimo
g = x-> -(-(0.4/((x^2+1))^(1/2))+((x^2+1)^(1/2))*(1-0.4/(x^2+1))-x) #funcao original

#definicao da primeira derivada
df(x) = 1.2*x/((x^2+1)^(3/2))+(1-0.4/(x^2+1))*x/((x^2+1)^(1/2))-1

#dois pontos iniciais
x0 = 3.4
x1 = 3.8

#convergencia
c_tol = 1.0e-8

as,bs = bracket_sign_change(df,x0,x1)
as,df(as),bs,df(bs)

af,bf = bisection(df,as,bs,c_tol)


println("x* = ", af)
println("df = ", -df(af))
println("f = ", g(af))

plot(g,as-1.,bs+1.,lw=2,draw_arrow="true",label="Function")
scatter!([as,bs],[g(as),g(bs)],label="interval")
scatter!([af,bf],[g(af),g(bf)],label="x*")
