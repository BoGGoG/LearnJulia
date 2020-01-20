# ToDo:
# - Box counting dimension
# - arbitray dimensions

using Random
using Plots
using BenchmarkTools
using ProgressBars

const numWalkers = 50000 # number of random walks
const N = 1000 # max number of steps in a single walk
const startPos = [0,0]
##

# task:
# 2D random walk, 4 directions (up, down, right, left)
# - find average distance R from origin after N steps
# - assume R ~ N^a, calc a


function randomStep(startPos)
    # from startPos perform one random step in direction N, S, E or W
    # return: position after this one step
    pos = copy(startPos)
    directions = ['N', 'S', 'E', 'W']
    dir = rand(directions)
    if dir == 'N'
        pos[2] += 1
    elseif dir == 'S'
        pos[2] -= 1
    elseif dir == 'W'
        pos[1] -= 1
    elseif dir == 'E'
        pos[1] += 1
    else
        print("Something went wron in randomStep")
    end
    return pos
end

function createRandomWalkHistory(startPos = [0,0], Nsteps = 10000)
    # Let random walker start from startPos and do Nsteps steps.
    # return: history of positions in form of array of size (2, Nsteps+1)
    # posArray[x/y, step]
    posArray = zeros(Int64, 2, Nsteps + 1)
    posArray[:, 1] = copy(startPos)
    for step in 1:1:Nsteps
        posArray[:, step + 1] = randomStep(posArray[:,step])
    end
    return posArray
end

function randomWalkEndPos(; startPos = [0,0], Nsteps = 10000)
    pos = copy(startPos)
    for step in 1:1:Nsteps
        pos = randomStep(pos)
    end
    return pos
end


function randomWalkers(; numWalkers = 100, Nsteps = 10000)
    endPosArray = zeros(Int64, 2, numWalkers)
    for walker in eachindex(endPosArray[1,:])
        endPosArray[:, walker] = randomWalkEndPos(Nsteps = Nsteps)
    end
    return endPosArray
end

function randomWalkers2(; numWalkers = 100, Nsteps = 10000)
    endPosArray = zeros(Int64, 2, numWalkers)
    endPosArray = mapslices(x -> randomWalkEndPos(Nsteps = Nsteps), endPosArray, dims = [1])
    return endPosArray
end


##

# plot one random walk
history = createRandomWalkHistory([1,1], N)
p = plot(history[1,:], history[2,:])

##

# plot end positions of K random walks
# endPosArray = randomWalkers(numWalkers = numWalkers, Nsteps = N)
# scatter(endPosArray[1,:], endPosArray[2,:])
# scatter(endPosArray2[1,:], endPosArray2[2,:])

########################################
# Statistical Evaluation
########################################
using Statistics

function radialDistance(x, y)
    return sqrt(x^2 + y^2)
end

endPosArray = randomWalkers2(numWalkers = numWalkers, Nsteps = N)
distancesArray = mapslices(vector -> radialDistance(vector[1], vector[2]), endPosArray, dims = [1])
histogram(distancesArray[:], normalize = false)

expecVal = mean(distancesArray)

function calcMeanDist(N; numWalkers = 10000)
    endPosArray = randomWalkers2(numWalkers = numWalkers, Nsteps = N)
    distancesArray = mapslices(vector -> radialDistance(vector[1], vector[2]), endPosArray, dims = [1])
    expecVal = mean(distancesArray)
    return expecVal
end

Nspan = 1:200
meanDistArr = calcMeanDist.(Nspan, numWalkers = 10000)

using LsqFit
@. model(N, p) = p[1]*N^p[2]

xData = Nspan
yData = meanDistArr
p0 = [0.5, 0.5]
fit = curve_fit(model, xData, yData, p0)

plot(meanDistArr, label = "'data'",
    xlabel = "N",
    ylabel = "E [dx]",
    legend = :bottomright)
plot!(model(xData, fit.param),
    label = "fit $(round(fit.param[1], digits = 2)) N^$(round(fit.param[2], digits = 2))")

println("The distance scales like N^$(round(fit.param[2], digits = 2))")
