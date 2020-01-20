using LinearAlgebra, Plots, ApproxFun

function chebPol(n::Int; spaceRange::Interval = -1..1)
    coefs = zeros(Int, n+1)
    coefs[n+1] = 1
    return Fun(Chebyshev(spaceRange), coefs)
end

##
# p = plot()
# spaceRange = 0..1
# plotPoints = spaceRange.left:0.05:spaceRange.right
# chebPol(4).(plotPoints)
# for n in 0:4
#     plot!(p, plotPoints, chebPol(n, spaceRange = 0..1).(plotPoints))
# end
# p
##
