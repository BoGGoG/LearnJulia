using DifferentialEquations
using Plots; gr()
using Benchmarktools

function henonHeiles!(du, u, p, t)
    # u[1] = x
    # u[2] = px
    # u[3] = y
    # u[4] = py
    du[1] = u[2]
    du[2] = -u[1] - 2*p[1]*u[1]*u[3]
    du[3] = u[4]
    du[4] = -u[3] - p[1]*(u[1]^2 - u[3]^2)
end

u0 = [0.5; 0.0; -0.5; 0.0]
tspan = (0.0,100.0)
p = [0.1]
prob = ODEProblem(henonHeiles!, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, vars = (1,3))
plot(sol)
