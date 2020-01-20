using DifferentialEquations

f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
using Plots
plot(
    sol,
    linewidth = 5,
    title = "Solution to the linear ODE with a thick line",
    xaxis = "Time (t)",
    yaxis = "u(t) (in μm)",
    label = "My Thick Line!",
) # legend=false
plot!(
    sol.t,
    t -> 0.5 * exp(1.01t),
    lw = 3,
    ls = :dash,
    label = "True Solution!",
)

gui()

# System of Equations
using ParameterizedFunctions

g! = @ode_def SystemOfEquations begin
    dx = σ * (y - x)
    dy = x * (ρ - z) - y
    dz = x * y - β * z
end σ ρ β

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
p = [10.0, 28.0, 8 / 3]
prob = ODEProblem(g!, u0, tspan, p)
sol = solve(prob)
plot(sol)
plot(sol, vars = (1, 2, 3))

# Nonhomogeneous ODEs
# ====================================
using DifferentialEquations
using Plots

l = 1.0
m = 1.0
g = 9.81

# pendulum! = @ode_def Pendulum begin
#     dθ = ω
#     dω = -3/2 * g/l * sin(θ) + 3/(m*l^2) * p(t)
function pendulum!(du, u, p, t)
    du[1] = u[2]
    du[2] = -3g / (2l) * sin(u[1]) + 3 / (m * l^2) * p(t)
end

θ0 = 0.01
ω0 = 0.0
u0 = [θ0, ω0]
tspan = (0.0, 10.0)

# M = t->0.1sin(t)
MM(t) = 0.1 * sin(t)
prob = ODEProblem(pendulum!, u0, tspan, MM)
sol = solve(prob)
gr()
plot(
    sol,
    linewidth = 2,
    xaxis = "t",
    label = ["θ [rad]" "ω [rad/s]"],
    layout = (2, 1),
)

# MATRIX equation, let u be a matrix
# =====================================
A = [
    1.0 0 0 -5
    4 -2 4 -3
    -4 0 0 1
    5 -2 2 3
]
u0 = rand(4, 2)
tspan = (0.0, 1.0)
function fMatrix(u, p, t)
    return A * u
end

prob = ODEProblem(f, u0, tspan)
sol = solve(prob)
plot(sol)

function fBetter!(du, u, p, t)
    mul!(du, A, u)
end

p = 0
prob = ODEProblem(fBetter!, u0, tspan, p)
sol = solve(prob)
plot(sol)

using LinearAlgebra
using StaticArrays

A = @SMatrix [
    1.0 0 0 -5
    4 -2 4 -3
    -4 0 0 1
    5 -2 2 3
]
u0 = @SMatrix rand(4, 2)
tspan = (0.0, 1.0)
prob = ODEProblem(fMatrix, u0, tspan)
sol = solve(prob)
plot(sol)
