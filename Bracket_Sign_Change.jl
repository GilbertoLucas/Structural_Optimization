


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
    

df(x) = 10*cos(10*x)-3*sin(3*x)

as,bs = bracket_sign_change(df,3.4,3.8)
as,df(as),bs,df(bs)