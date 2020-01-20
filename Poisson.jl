# Poisson Equation
# =======================
using ParameterizedFunctions
f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
gD(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)

f  = @fem_def((x),TestF,begin
  sin(α.*x).*cos(α.*y)
end,α=>6.28)
gD = @fem_def (x) TestgD begin
  sin(α.*x).*cos(α.*y)/β
end α=>6.28) β=>79.0
