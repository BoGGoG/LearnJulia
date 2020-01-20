include("myRandomWalk.jl")
using Plots
import .myRandomWalk

const numWalkers = 1000
const N = 100 # max number of steps in a single walk
const startPos = [0,0]

myRandomWalk.say_hello()
myRandomWalk.randomStep([0,0])

history = myRandomWalk.createRandomWalkHistory([1,1], N)
p = myRandomWalk.plot(history[1,:], history[2,:])
p

endPosArray = myRandomWalk.randomWalkers(numWalkers = numWalkers, Nsteps = N)
distancesArray = mapslices(vector -> myRandomWalk.radialDistance(vector[1], vector[2]), endPosArray, dims = [1])
myRandomWalk.histogram(distancesArray[:], normalize = false)


const Nspan = 1:200
print(typeof(Nspan))
meanDistArr = myRandomWalk.calcMeanDist.(Nspan, numWalkers = 10000)

fitResults = myRandomWalk.fitMeanDist(1:30; numWalkers = 1000)
myRandomWalk.plotFitMeanDist(fitResults)

myRandomWalk.randomWalkUntilBorder([0,0], 5)
myRandomWalk.randomCityWalkers(10, 20)

myRandomWalk.calcCityMeanSteps(20)
myRandomWalk.calcAndPlotCityMeanSteps(1:20; numberOfWalkers = 1000, randomStart = true)
myRandomWalk.calcAndPlotCityMeanSteps!(1:20; numberOfWalkers = 1000, randomStart = false)
