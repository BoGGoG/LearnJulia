using DifferentialEquations
using ParameterizedFunctions
function f(du, u, p, t)
    x, y = u
    α, β, γ, δ = p
    du[1] = α*x - β*x*y
    du[2] = -γ*y + δ*x*y
end

p = (1.5, 1.0, 3.0, 1.0); u0 = [1.0; 1.0]
tspan = (0.0, 10.0)
prob1 = ODEProblem(f, u0, tspan, p)
sol = solve(prob1)
using Plots
plot(sol)

g = @ode_def LotkaVolterra begin
    dx = α*x - β*x*y
    dy = -γ*y + δ*x*y
end α β γ δ


prob2 = ODEProblem(g, u0, tspan, p)
sol2 = solve(prob2)
plot(sol2)


function noise(du, u, p, t)
    du[1] = 0.2u[1]
    du[2] = 0.2u[2]
end
prob3 = SDEProblem(g, noise, u0, tspan, p)
sol3 = solve(prob2)
plot(sol3)
