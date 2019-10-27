#include <iostream>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

//Some parameters for the simulation
const int POP_SIZE = 100; //population size
const int SIZE = 40; //The number of members in each tree
const int MAX_DEGREE = 4; //The degree constraint
const int MAX_EDGE_WEIGHT = 10; //Maximum weight for a single edge
const int MUTATION_PROBABILITY = 100; //Reciprocal of probability for mutation.
const int NUM_GENERATIONS = 1000;
const int RUNS_PER_INSTANCE = 30;
const int NUM_INSTANCES = 10;

//Function declarations
void pruferString(int[],int, int);
void displayArray(int[],int);
int getFitness(int[],int, int[][SIZE]);
void blobDecode(int[], int, int[][2]);
bool path(int,int,int[][2],int[]);
void displayArray(double arr[], int size);
void generateNeighbors(int [SIZE-2], int [2*(SIZE-2)][SIZE-2], int );
bool meetsConstraint(int [], int);


int main()
{
	//Seeding for random number generation
	srand((unsigned int)time(NULL));

	//Data collection arrays
	double averages[NUM_GENERATIONS];
	int minima[NUM_GENERATIONS];
	double minAverages[RUNS_PER_INSTANCE];
	int globalMinima[RUNS_PER_INSTANCE];
	double instanceMinAverages[NUM_INSTANCES];
	int instanceGlobalMinima[NUM_INSTANCES];
	int tabuMinima[NUM_GENERATIONS];
	int tabuRunMinima[RUNS_PER_INSTANCE];
	int instanceTabuMinima[NUM_INSTANCES];

	ofstream output1("averages.dat");
	ofstream output2("minima.dat");
	ofstream output3("minAverages.dat");
	ofstream output4("globalMinima.dat");
	ofstream output5("instanceMinAverages.dat");
	ofstream output6("instanceGlobalMinima.dat");
	ofstream output7("runTabuMinima.dat");
	ofstream output8("instanceTabuMinima.dat");
	
	//Create a random adjacency list for the connected graph
	int adjacencyList[SIZE][SIZE];
	int instanceNum = 0;
	while (instanceNum < NUM_INSTANCES)
	{
		for(int i = 0; i < SIZE; i++)
		{
			adjacencyList[i][i] = 0; //no vertex is connected to itself
			for(int j = i+1; j < SIZE; j++)
			{
				adjacencyList[i][j] = rand() % MAX_EDGE_WEIGHT + 1; //+1 ensures weight > 0 for i != j, so the graph is connected
				adjacencyList[j][i] = adjacencyList[i][j]; //non-directional adjacency lists are always symmetric
			}
		}
		int runNum = 0;
		while(runNum < RUNS_PER_INSTANCE)
		{
			//Initialize a population of Prufer strings that meet the 
			//degree constraint requirements (set in MAX_DEGREE constant)
			int population[POP_SIZE][SIZE - 2];	
			for (int i = 0; i < POP_SIZE; i++)
			{
				pruferString(population[i], SIZE, MAX_DEGREE);
			}
			
			//Set initial Tabu Prufer string to first element of generated population.
			int currentTabu[SIZE-2];
			int tabuList[SIZE-2];
			for (int i = 0; i < SIZE-2; i++)
			{
				currentTabu[i] = population[0][i];
				tabuList[i] = 0;
			}
			//WHILE (GENERATION < TOTALGENERATIONS)
			int genNum = 0;
			while(genNum < NUM_GENERATIONS)
			{
				//Get some basic statistics about the population
				cout << "Instance Number: " << instanceNum + 1 << endl;
				cout << "Run Number: " << runNum + 1 << endl;
				cout << "Generation Number: " << genNum + 1 << endl;
				int sum = 0;
				double avg = 0;
				int tabuFitness = getFitness(currentTabu, SIZE, adjacencyList);
				int minElement = getFitness(population[0],SIZE,adjacencyList);
				for (int i = 0; i<POP_SIZE; i++)
				{
					int c = getFitness(population[i],SIZE,adjacencyList);
					sum += c;
					if (minElement > c)
						minElement = c;
				}
				avg = sum / POP_SIZE;
				averages[genNum] = avg;
				minima[genNum] = minElement;
				cout << "Average fitness: " << avg << endl;	
				cout << "Minimum: " << minElement << endl;
				cout << "Tabu Fitness: " << tabuFitness << endl;
				//1. CONDUCT TWO-TOURNAMENT SELECTION OF PARENTS
				int parentPopulation[POP_SIZE][SIZE-2];
				for (int i = 0; i < POP_SIZE; i++)
				{
					//1a. SELECT TWO PARENTS AT RANDOM
					int parentIndex1 = rand() % POP_SIZE;
					int parentIndex2 = rand() % POP_SIZE;
					//1b. PLACE PARENT WITH HIGHER FITNESS INTO PARENT POPULATION
					if (getFitness(population[parentIndex1],SIZE,adjacencyList)<=getFitness(population[parentIndex2],SIZE,adjacencyList))
					{	
						for(int j = 0; j< SIZE - 2; j++)
							parentPopulation[i][j] = population[parentIndex1][j];
					}
					else
					{
						for(int j = 0; j < SIZE - 2; j++)
							parentPopulation[i][j] = population[parentIndex2][j];
					}
				}

				//2. RECOMBINE PARENTS TO GET CHILD POPULATION
				int childPopulation[POP_SIZE][SIZE-2];	
				for (int i = 0; i < POP_SIZE; i++)
				{
					//2a.SELECT TWO RANDOM PARENTS FROM THE PARENT POPULATION
					int parentIndex3 = rand() % POP_SIZE;
					int parentIndex4 = rand() % POP_SIZE;
					//2b.AUGMENT PARENTAL CHROMOSOMES WITH DEGREE COUNT
					int count1[SIZE];
					int count2[SIZE];
					int maxCount[SIZE];
					for (int k = 0; k<SIZE; k++)
					{
						int c1 = 0;
						int c2 = 0;
						for (int j = 0; j<SIZE-2;j++)
						{
							if (k==parentPopulation[parentIndex3][j])
								c1++;
							if (k==parentPopulation[parentIndex4][j])
								c2++;
						}
						count1[k] = c1;
						count2[k] = c2;
						maxCount[k] = max(c1,c2);
					}
					//2c.CONDUCT MODIFIED CYCLE CROSSOVER 
					int randSpot1 = rand() % (SIZE-2);
					int randSpot2 = rand() % (SIZE-2);
					if (randSpot1 == randSpot2)
						randSpot2 = (randSpot2 + 1) % (SIZE-2);
					int c1 = min(randSpot1,randSpot2);
					int c2 = max(randSpot1,randSpot2);
					int newChromosome[SIZE-2];
					for (int k = c1; k <= c2; k++)
					{
						newChromosome[k]= parentPopulation[parentIndex3][k];
						maxCount[newChromosome[k]]--;
					}
					int j = ((c2+1)%(SIZE-2));
					int currentParentIndex = (c2+1)%(SIZE-2);
					while (j%(SIZE-2) != c1)
					{
						if (maxCount[parentPopulation[parentIndex4][currentParentIndex % (SIZE-2)]] > 0)
						{
							newChromosome[j%(SIZE-2)] = parentPopulation[parentIndex4][currentParentIndex % (SIZE-2)];
							j++;
							currentParentIndex++;
						}
						else
							currentParentIndex++;
					}
					//2d.PLACE CHILD CHROMOSOME INTO CHILD POPULATION
					for (int j = 0; j < SIZE-2;j++)
					{
						childPopulation[i][j] = newChromosome[j];
					}
				}
				//3. MUTATE CHILD POPULATION
				for (int i = 0; i < POP_SIZE; i++)
				{
					for(int j = 0; j < SIZE-2; j++)
					{
						int randInt = rand() % MUTATION_PROBABILITY; 
						if (randInt == 1)
						{	
							int count = 0;
							for(int k = 0; k < SIZE-2; k++)
							{
								if ((childPopulation[i][j]+1)%SIZE==childPopulation[i][k])
									count++;
							}
							if (count < MAX_DEGREE)
							{
								childPopulation[i][j] = (childPopulation[i][j]+1)%SIZE;
							}
						}
					}
				}
				//4. SET POPULATION TO CHILD POPULATION
				for (int i = 0; i < POP_SIZE; i++)
				{
					for (int j = 0; j < SIZE-2; j++)
					{
						population[i][j] = childPopulation[i][j];
					}
				}
				
				//5. RUN TABU SEARCH ITERATION
				//Generate neighbor list of current Tabu
				int neighborList[2*(SIZE-2)][SIZE-2];
				generateNeighbors(currentTabu, neighborList, SIZE);
				//Replace currentTabu with the best neighbor that meets constraint requirements
				int replaceIndex = 0;
				for (int z = 0; z < 2*(SIZE-2); z++)
				{
					if((getFitness(neighborList[replaceIndex], SIZE, adjacencyList)>getFitness(neighborList[z],SIZE, adjacencyList)) && (tabuList[z%(SIZE-2)]<=0) && meetsConstraint(neighborList[replaceIndex], SIZE))
					{
						replaceIndex = z;
					}
				}
				tabuList[replaceIndex%(SIZE-2)] = 35;
				for (int x = 0; x < SIZE-2; x++)
					tabuList[x]--;
				for (int z = 0; z < SIZE-2; z++)
				{
					currentTabu[z] = neighborList[replaceIndex][z];
				}
				displayArray(currentTabu, SIZE-2);
				tabuMinima[genNum]=getFitness(currentTabu, SIZE, adjacencyList);
				//Add operation to the Tabu list
				//6. GENERATION += 1, REPEAT WHILE LOOP
				genNum++;
			}
			double minimumAverage = averages[0];
			int globalMinimum = minima[0];
			int tabuRunMinimum = tabuMinima[0];	
			for(int i = 0; i < NUM_GENERATIONS; i++)
			{
				if (minimumAverage > averages[i])
					minimumAverage = averages[i];
				if (globalMinimum > minima[i])
					globalMinimum = minima[i];	
				if (tabuRunMinimum > tabuMinima[i])
					tabuRunMinimum = tabuMinima[i];
				output1 << averages[i] << "\n";
				output2 << minima[i] << "\n";
			}
			minAverages[runNum] = minimumAverage;
			globalMinima[runNum] = globalMinimum;
			tabuRunMinima[runNum] = tabuRunMinimum;
			output3 << minAverages[runNum]<< endl;
			output4 << globalMinima[runNum] << endl;
			output7 << tabuRunMinima[runNum] << endl;
			runNum++;
		}
		double instanceMinAverage = minAverages[0];
		int instanceGlobalMinimum = globalMinima[0];
		int tabuInstanceMinimum = tabuRunMinima[0];
		for (int i = 0; i < RUNS_PER_INSTANCE; i++)
		{
			if (instanceMinAverage > minAverages[i])
				instanceMinAverage = minAverages[i];
			if (instanceGlobalMinimum > globalMinima[i])
				instanceGlobalMinimum = globalMinima[i];
			if (tabuInstanceMinimum > tabuRunMinima[i])
				tabuInstanceMinimum = tabuRunMinima[i];
		}
		output5 << instanceMinAverage << endl;
		output6 << instanceGlobalMinimum << endl;
		output8 << tabuInstanceMinimum << endl;
		instanceNum++;
	}
	output1.close();
	output2.close();
	output3.close();
	output4.close();
	output5.close();
	output6.close();
	output7.close();

	return 0;
}

void pruferString(int arr[],int size, int degree)
{
	int degreeTracker[size];
	for(int i = 0; i < size; i++)
	{
		degreeTracker[i] = degree;
	}
	
	int j = 0;	
	while(j < size - 2)
	{
		int randomInt = rand() % size;
		if(degreeTracker[randomInt] > 0)
		{
			arr[j] = randomInt;
			degreeTracker[randomInt] -= 1;
			j++;
		}	
	}

}

bool meetsConstraint(int pruferString[], int size)
{
	for (int i = 0; i < size; i++)
	{
		int currentDegree = 0;
		for(int j = 0;j<size-2;j++)
		{
			if (pruferString[j] == i)
				currentDegree++;
		}
		if (currentDegree > MAX_DEGREE)
			return false;
	}
	return true;
}

void displayArray(int arr[], int size)
{
	for(int i = 0; i < size; i++)
	{
		cout << arr[i] << " ";
	}
	cout << endl;
}

void displayArray(double arr[], int size)
{
	for(int i = 0; i < size; i++)
	{
		cout << arr[i] << " ";
	}
	cout << endl;
}

int getFitness(int pruferString[], int size, int adjacencyMatrix[][SIZE])
{
	int fitness = 0;
	//GET EDGE LIST FROM PRUFER STRING INPUT USING BLOB DECODE
	int edgeList[size - 1][2];
	for (int i = 0; i < size -1; i++)
	{
		edgeList[i][0] = -1;
		edgeList[i][1] = -1;
	}
	blobDecode(pruferString, size, edgeList);
	//SUM EDGE LIST VALUES FROM ADJACENCY MATRIX
	for (int j = 0; j < size -2; j++)
	{
		fitness += adjacencyMatrix[edgeList[j][0]][edgeList[j][1]];
	}
	return fitness;
}

void blobDecode(int pruferString[], int size, int edgeList[][2])
{
	int blob[size];
	blob[0] = 0;
	for(int i=1;i<=size-1;i++)
		blob[i]=1;
	int blobPointer = 0;
	for(int i=0;i<=size-2;i++)
	{
		blob[i+1]=0;
		if (path(pruferString[i], size, edgeList, blob))
		{
			edgeList[i][0] = i+1;
			edgeList[i][1] = pruferString[i];
		}
		else
		{
			edgeList[i][0] = i+1;
			edgeList[i][1] = blobPointer;
			blobPointer = pruferString[i];
		}
	}
	edgeList[size-1][0]=size-1;
	edgeList[size-1][1]=blobPointer;
}

bool path(int a,int size,int edgeList[][2],int blob[])
{
	bool atEnd = false;
	int currentVertex = a;
	int previousVertex = a;
	while(!atEnd)
	{
		previousVertex = currentVertex;
		if (blob[currentVertex] == 1)
			return true;
		for (int i = 0; i < size - 2; i++)
		{
			if(edgeList[i][0]==currentVertex)
			{
				currentVertex = edgeList[i][1];
				break;
			}
		}
		if(currentVertex == previousVertex)
			atEnd = true;
	}
	return false;
}

void generateNeighbors(int pruferString[SIZE-2], int neighborList[2*(SIZE-2)][SIZE-2], int size)
{
	int currentNeighbor[size-2];
	for (int i = 0; i < 2*(size-2); i++)
	{
		for (int j = 0; j < size-2;j++)	
		{
			if((j%(size-2)) != (i%(size-2)))
				currentNeighbor[j] = pruferString[j];
			else
			{
				if (i < size-2)
					currentNeighbor[j] = (pruferString[j]+1)%(size);
				else
					currentNeighbor[j] = (pruferString[j]-1 + size)%(size);
			}

		}
		for (int k = 0; k<size-2;k++)
		{
			neighborList[i][k] = currentNeighbor[k];
		}
	}
}

