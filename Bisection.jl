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
    

    
f = x->sin(10x)+cos(3x)
df(x) = 10*cos(10*x)-3*sin(3*x)

as,bs = bracket_sign_change(df,3.4,3.8)
as,df(as),bs,df(bs)

af,bf = bisection(df,as,bs,1.0e-8)

df(af),f(af)