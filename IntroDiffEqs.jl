# Following https://tutorials.juliadiffeq.org/html/introduction/01-ode_introduction.html
using DifferentialEquations
using Plots; gr()
f(u,p,t) = 0.98u
u0 = 1.0
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob)
gr()
plot(sol, linewidth = 1, title = "u(t)'=0.98u(t)",
    xaxis = "Time (t)", yaxis="u(t)", label = "Numerical Solution")
plot!(sol.t, t->1.0*exp(0.98t), lw=2, ls=:dash, label = "True Solution")

# better tollerance
sol2 = solve(prob, abstol = 1e-8, reltol = 1e-8)
sol2.t
plot!(sol2, linewidth = 1, label = "Numerical Solution, abstol 1e-8")

# Systems of ODEs: The Lorenz Equation
function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ*(u[2] - u[1])
    du[2] = u[1]*(ρ - u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end
u0 = [1,0, 0.0, 0.0]
p = (10, 28, 8/3)
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob)

plot(sol)
plot(sol, vars = (1,2,3))

# Lotka Volterra using ParameterizedFunctions
function lotka_volterra!(du, u, p, t)
    du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
    du[2] = -p[3]*u[2] + p[4]*u[1]*u[2]
end

using ParameterizedFunctions
lv! = @ode_def LotkaVolterra begin
    dx = a*x - b*x*y
    dy = -c*y + d*x*y
end a b c d

u0 = [1.0, 1.0]
p = (1.5, 1.0, 3.0, 1.0)
tspan = (0.0, 10.0)
prob2 = ODEProblem(lv!, u0, tspan, p)
sol2 = solve(prob2)
plot(sol2)

print(lv!.Jex)

using Latexify
print(latexify(lv!))
