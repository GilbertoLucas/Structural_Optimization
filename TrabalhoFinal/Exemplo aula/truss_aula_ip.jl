import JuMP as J
using Ipopt

# Dados de entrada 

Sy=350.0;
r1=0.4;
E=2000.0;
P=200.0;
B=1.0;


f(r2,H) = pi*(r2^2)*(H^2+B^2)^(1/2)

g1(r2,H) = P*B/(pi*H*(r1^2)*Sy) - 1;
g2(r2,H) = P*((B^2+H^2)^(1/2))/(H*pi*r2^2*Sy)-1;
g3(r2,H) = 4*P*(B^2+H^2)^(3/2)/(H*pi^3*r2^4*E)-1;



model = J.Model(Ipopt.Optimizer)

J.set_optimizer_attribute(model, "print_level",0)

#registrando o número de variáveis 

J.register(model, :g1, 2, g1; autodiff = true)
J.register(model, :g2, 2, g2; autodiff = true)
J.register(model, :g3, 2, g3; autodiff = true)
J.register(model, :f, 2, f; autodiff = true)

J.@variable(model, 0.4 <= r2 <= 0.8);
J.@variable(model, 1.0 <= H <= 6.0);


J.@NLobjective(model, Min, f(r2,H));

J.@NLconstraint(model, g1(r2,H) <= 0);
J.@NLconstraint(model, g2(r2,H) <= 0);
J.@NLconstraint(model, g3(r2,H) <= 0);


println(model)
J.optimize!(model)

println("Termination status: ", J.termination_status(model))
println("Primal status: ", J.primal_status(model))
println("f*: ", J.objective_value(model))
println("r2* : ", J.value(r2), ", H* : ", J.value(H))


a = 4;