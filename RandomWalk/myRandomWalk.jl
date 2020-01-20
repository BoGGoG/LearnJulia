module myRandomWalk
export say_hello, randomStep, createRandomWalkHistory, randomWalkers, radialDistance, calcMeanDist, fitMeanDist, plotFitMeanDist, randomWalkUntilBorder, randomCityWalkers, calcCityMeanSteps, calcAndPlotCityMeanSteps

using Random
using Plots
using Statistics
using LsqFit

say_hello() = println("hello")
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
        print("Something went wrong in randomStep")
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
    endPosArray = mapslices(x -> randomWalkEndPos(Nsteps = Nsteps), endPosArray, dims = [1])
    return endPosArray
end

function randomWalkersForLoop(; numWalkers = 100, Nsteps = 10000)
    endPosArray = zeros(Int64, 2, numWalkers)
    for walker in eachindex(endPosArray[1,:])
        endPosArray[:, walker] = randomWalkEndPos(Nsteps = Nsteps)
    end
    return endPosArray
end

function radialDistance(x, y)
    return sqrt(x^2 + y^2)
end

function calcMeanDist(N; numWalkers = 10000)
    endPosArray = randomWalkers(numWalkers = numWalkers, Nsteps = N)
    distancesArray = mapslices(vector -> radialDistance(vector[1], vector[2]), endPosArray, dims = [1])
    expecVal = mean(distancesArray)
    return expecVal
end

"""Calculate the mean walked distance for Nspan steps
    by taking the mean of numWalkers walkers for every
    number of steps and fitting to N^a

    return: fit, (Nspan, meanDistArr)"""
function fitMeanDist(Nspan::UnitRange; numWalkers = 10000)
    @. model(N, p) = p[1]*N^p[2]

    meanDistArr = calcMeanDist.(Nspan, numWalkers = 10000)
    xData = Nspan
    yData = meanDistArr
    p0 = [0.5, 0.5]
    fit = curve_fit(model, xData, yData, p0)
    return fit, (Nspan, meanDistArr)
end

function plotFitMeanDist(fitResults)
    fit, tup = fitResults
    (Nspan, meanDistArr) = tup

    @. model(N, p) = p[1]*N^p[2]
    xData = Nspan
    yData = meanDistArr

    plot(meanDistArr, label = "'data'",
        xlabel = "N",
        ylabel = "E [dx]",
        legend = :bottomright)
    plot!(model(xData, fit.param),
        label = "fit $(round(fit.param[1], digits = 2)) N^$(round(fit.param[2], digits = 2))")
end

"""if outside of his box then he's done
numberOfSites: square lattice, so in each direction,
e.g. numberOfSites = 4 --> -2, -1, 0, 1, 2
e.g. numberOfSites = 5 --> -2, -1, 0, 1, 2
e.g. numberOfSites = 6 --> -3, -2, -1, 0, 1, 2, 3
"""
function isOutside(pos::Array{Int64,1}, numberOfSites::Int64)
    maxPos = [numberOfSites÷2, numberOfSites÷2]
    return any(abs.(pos).> maxPos)
end

function randomWalkUntilBorder(startPos::Array{Int64,1}, numberOfSites::Int64; maxSteps = 1000000)
    stepsUntilOutside = 0
    pos = copy(startPos)
    for step in 1:maxSteps
        pos = randomStep(pos)
        if isOutside(pos, numberOfSites)
            break
        end
        stepsUntilOutside = step
    end
    return stepsUntilOutside, pos
end

function randomPos(numberOfSites::Int64)
    maxPos = [numberOfSites÷2, numberOfSites÷2]
    x = rand(-maxPos[1]:maxPos[1])
    y = rand(-maxPos[2]:maxPos[2])
    return [x,y]
end

"""let walkers walk in square city of floor(numberOfSites) sites
until boudary is reached.
Return: Number of steps"""
function randomCityWalkers(numberOfWalkers::Int64, numberOfSites::Int64; randomStart = true)
    stepsArray = Array{Int64}(undef, numberOfWalkers)
    if randomStart
        map!(x -> randomWalkUntilBorder(randomPos(numberOfSites), numberOfSites)[1], stepsArray, stepsArray)
    else
        map!(x -> randomWalkUntilBorder([0,0], numberOfSites)[1], stepsArray, stepsArray)
    end
    return stepsArray
end

function calcCityMeanSteps(numberOfSites::Int; numberOfWalkers = 10000, randomStart = true)
    return mean(randomCityWalkers(numberOfWalkers, numberOfSites, randomStart = randomStart))
end

function calcAndPlotCityMeanSteps(Nspan::UnitRange; numberOfWalkers = 10000, randomStart = true)
    meanDistArray = calcCityMeanSteps.(Nspan; numberOfWalkers = numberOfWalkers, randomStart = randomStart)
    plot(meanDistArray)
end

function calcAndPlotCityMeanSteps!(Nspan::UnitRange; numberOfWalkers = 10000, randomStart = true)
    meanDistArray = calcCityMeanSteps.(Nspan; numberOfWalkers = numberOfWalkers, randomStart = randomStart)
    plot!(meanDistArray)
end

end
