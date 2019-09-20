import random
import numpy as np
import matplotlib.pyplot as plt
import pylab

def main():
    #Basic population variables
    populationSize = 100
    generationNumber = 200
    mutationProbability = 100

    #Lists for plotting
    averageFitnesses = []
    maximumFitnesses = []

    #Random instance generation
    matrixSize = random.randint(30,50)
    instanceMatrix = np.zeros((matrixSize, matrixSize))
    for i in range(0,matrixSize):
        for j in range(0, matrixSize):
            instanceMatrix[i][j] = random.randint(-100,100)
    print("Instance Matrix:")
    print(instanceMatrix)

    #Initialize population
    population = initializePopulation(populationSize, matrixSize)
    print("Initial population:")
    #displayGenePool(population, instanceMatrix)

    #Selection, recombination, and mutation loop
    currentGeneration = 1
    while currentGeneration <= generationNumber:

        #Obtain parent pool from current population using 2-tournament selection
        parentPool = getParentPool(population, instanceMatrix, populationSize)
        #print("Parent pool:")
        #displayGenePool(parentPool, instanceMatrix)

        #Randomly recombine pairs of parents to get new pool
        population = getRecombinedPool(parentPool, populationSize, matrixSize)
        #print("Non-mutated population: ")
        #displayGenePool(population, instanceMatrix)

        #Randomly mutate the members of the new population
        population = getMutatedPool(population, mutationProbability)
        #print("Mutated population: ")
        #displayGenePool(population, instanceMatrix)

        #Display some fitness statistics
        print("Generation: %d" % currentGeneration)
        displayFitnessStats(population, instanceMatrix)

        #Update statistic lists for plotting    
        fitnessList = []
        for i in population:
            fitnessList.append(getFitness(i, instanceMatrix))
        averageFitnesses.append(sum(fitnessList)/len(fitnessList))
        maximumFitnesses.append(max(fitnessList))

        #Go to next generation
        currentGeneration += 1
        
    plt.plot(range(1, generationNumber+1), maximumFitnesses, 'ro',label = "Maximum Fitness")
    plt.plot(range(1, generationNumber+1) ,averageFitnesses, 'bo', label = "Average Fitness")
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    pylab.legend(loc='lower right')
    plt.show()

def displayGenePool(pool, instanceMatrix):
    for i in pool:
        print(i, end="\t\t")
        print(getFitness(i, instanceMatrix))
    print()
        

def generateBitList(matrixSize):
    bitList = []
    for i in range(0, matrixSize):
        randomBit = random.randint(0,1)
        bitList.append(randomBit)
    return bitList
    

def initializePopulation(populationSize, matrixSize):
    population = []
    for i in range(0, populationSize):
        population.append(generateBitList(matrixSize))
    return population

def getFitness(bitList, instanceMatrix):
    fitness = 0
    for i in range(0, len(bitList)):
        for j in range(0, len(bitList)):
            fitness += instanceMatrix[i][j]*bitList[i]*bitList[j]
    return fitness

def getParentPool(population, instanceMatrix, populationSize):
    parentPool = []
    for i in range(0,populationSize):
        randomParents = random.sample(population, 2)
        if getFitness(randomParents[0], instanceMatrix) >= getFitness(randomParents[1],instanceMatrix):
            parentPool.append(randomParents[0])
        else:
            parentPool.append(randomParents[1])
    return parentPool

def getRecombinedPool(parentPool, populationSize, matrixSize):
    childPool = []
    for i in range(0, populationSize):
        #Get two random parents
        randomIndices = random.sample(range(0, len(parentPool)-1), 2)
        parent1 = parentPool[randomIndices[0]]
        parent2 = parentPool[randomIndices[1]]

        #Recombine parent chromosomes at a random crossover point
        randomCrossoverPoint = random.randint(1, matrixSize)
        chrFirstHalf = parent1[0:randomCrossoverPoint]
        chrSecondHalf = parent2[randomCrossoverPoint:]
        child = chrFirstHalf
        for i in chrSecondHalf:
            child.append(i)
        childPool.append(child)
    return childPool

def getMutatedPool(pool, mutationProbability):
    for i in range(0, len(pool)):
        j = pool[i]
        for k in range(0, len(j)):
            if random.randint(1, mutationProbability) == 1:
                if j[k] == 0:
                    j[k] = 1
                else:
                    j[k] = 0
            pool[i] = j
    return pool

def displayFitnessStats(pool, instanceMatrix):
    fitnessList = []
    for i in pool:
        fitnessList.append(getFitness(i, instanceMatrix))
    print("Maximum fitness: %d" % max(fitnessList))
    averageFitness = sum(fitnessList)/len(fitnessList)
    print("Average fitness: %d" % averageFitness)


        
    
            
        

main()
