# Two dimensional projectile motion, no air resistance
# EoM:
# x(t) = x_0 + v0 cos(Θ)t
# y(t) = y_0 + v0 sin(Θ)t - g t^2 / 2
ENV["GKS_ENCODING"]="utf-8"
using Plots
gr()

g = 0.81
Θ0 = 45.0/180.0
v0 = 100.0
x0 = 0.0
y0 = 0.0
δt = 0.1
tspan = 0:100
times = tspan.start:δt:tspan.stop

function xyProjectile(x0, y0, v0, Θ0, t)
        x = x0 + v0*cos(Θ0)*t
        y = y0 + v0*sin(Θ0)*t - g*t^2 / 2.0
        return x,y
end

projectileMotion = xyProjectile.(x0, y0, v0, Θ0, times)

plot(projectileMotion,
        xlabel = "x [m]",
        ylabel = "y [m]",
        title = "Projectile Motion",
        label = "Θ0 = $Θ0, v0 = $v0")

using Roots
heightProjectile(t) = xyProjectile(x0, y0, v0, Θ0, t)[2]
landingTime = find_zeros(heightProjectile, 0, 1000)[2]
landingX = xyProjectile(x0, y0, v0, Θ0, landingTime)[1]

# calculating perfect angle

function calcRange(Θ)
        function calcy(y0, v0, Θ0, t)
                y = y0 + v0*sin(Θ0)*t - g*t^2 / 2.0
                return y
        end
        function calcx(x0, y0, v0, Θ0, t)
                x = x0 + v0*cos(Θ0)*t
                return x
        end
        groundTimes = find_zeros(t -> calcy(y0, v0, Θ, t), 0, 10000)
        return calcx(x0, y0, v0, Θ, groundTimes[end])
end


function calcΘopt(Θspan::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})
        return Θspan[findmax(calcRange.(Θspan))[2]]
end
calcΘopt(Θspan)

eps = 0.001
Θspan = 0:eps:π
Θopt = calcΘopt(Θspan)
xopt = calcRange(Θopt)
println("The optimal angle is $Θopt up to a precision of +-$eps")
plot(Θspan, calcRange.(Θspan))
scatter!([Θopt], [xopt], label = "($Θopt, $(round(xopt, digits = 2)))")

print([Θopt, calcRange(Θopt)])
