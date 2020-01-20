using Plots

function randomStep(startPos::Array{Int64,1})
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

##
# stepsArray = randomCityWalkers(10000, 15, randomStart = true)
# histogram(stepsArray)

##
##################################################
# Different city sizes
##################################################
using Statistics
using Distributed

function meanDist(numberOfSites::Int; numberOfWalkers = 10000)
    return mean(randomCityWalkers(numberOfWalkers, numberOfSites))
end

function meanDistOfSizes(sizes::Int; numberOfWalkers = 10000)
    return meanDist.(sizes, numberOfWalkers = numberOfWalkers)
end

function meanDistOfSizesParallel(sizes; numberOfWalkers = 10000)
    # return meanDist.(sizes, numberOfWalkers = numberOfWalkers)
    arr = collect(sizes)
    arr = pmap(n -> meanDist(n, numberOfWalkers = numberOfWalkers), arr)
end

meanDist(15, numberOfWalkers = 50000)
meanDistArray = meanDistOfSizes.(1:30, numberOfWalkers = 50000)
meanDistArray = meanDistOfSizesParallel(1:30, numberOfWalkers = 50000)
plot(meanDistArray)
