#include <iostream>
#include <string>
#include <time.h>
#include <stdlib.h>
using namespace std;

//Some parameters for the simulation
const int POP_SIZE = 100; //population size
const int SIZE = 10; //The number of members in each tree
const int MAX_DEGREE = 7; //The degree constraint
const int MAX_EDGE_WEIGHT = 10; //Maximum weight for a single edge
const int MUTATION_PROBABILITY = 100; //Reciprocal of probability for mutation.

//Function declarations
void pruferString(int[],int, int);
void displayArray(int[],int);
int getFitness(int[],int, int[][SIZE]);
void blobDecode(int[], int, int[][2]);
bool path(int,int,int[][2],int[]);

int main()
{
	//Seeding for random number generation
	srand((unsigned int)time(NULL));

	//Create a random adjacency list for the connected graph
	int adjacencyList[SIZE][SIZE];
	for(int i = 0; i < SIZE; i++)
	{
		adjacencyList[i][i] = 0; //no vertex is connected to itself
		for(int j = i+1; j < SIZE; j++)
		{
			adjacencyList[i][j] = rand() % MAX_EDGE_WEIGHT + 1; //+1 ensures weight > 0 for i != j, so the graph is connected
			adjacencyList[j][i] = adjacencyList[i][j]; //non-directional adjacency lists are always symmetric
		}
	}

	//Initialize a population of Prufer strings that meet the 
	//degree constraint requirements (set in MAX_DEGREE constant)
	int population[POP_SIZE][SIZE - 2];	
	for (int i = 0; i < POP_SIZE; i++)
	{
		pruferString(population[i], SIZE, MAX_DEGREE);
	}


	//WHILE (GENERATION < TOTALGENERATIONS)
	//1. CONDUCT TWO-TOURNAMENT SELECTION OF PARENTS
	int parentPopulation[POP_SIZE][SIZE-2];
	int sum = 0;
	for (int i = 0; i < POP_SIZE; i++)
	{
		//1a. SELECT TWO PARENTS AT RANDOM
		int parentIndex1 = rand() % POP_SIZE;
		int parentIndex2 = rand() % POP_SIZE;
		//1b. PLACE PARENT WITH HIGHER FITNESS INTO PARENT POPULATION
		if (getFitness(population[parentIndex1],SIZE,adjacencyList)>=getFitness(population[parentIndex2],SIZE,adjacencyList))
		{	
			for(int j = 0; j< SIZE - 2; j++)
				parentPopulation[i][j] = population[parentIndex1][j];
		}
		else
		{
			for(int j = 0; j < SIZE - 2; j++)
				parentPopulation[i][j] = population[parentIndex2][j];
		}

		for (int j = 0; j< POP_SIZE; j++)
		{
			sum += getFitness(population[j], SIZE, adjacencyList);
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
		if (i == POP_SIZE-1)
		{
			cout << "c1 is " << c1 << endl;
			cout << "c2 is " << c2 << endl;
			displayArray(parentPopulation[parentIndex3], SIZE-2);
			displayArray(parentPopulation[parentIndex4], SIZE-2);
		}
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
		if (i == POP_SIZE-1)
		{
			displayArray(newChromosome, SIZE-2);
		}	
		//2d.PLACE CHILD CHROMOSOME INTO CHILD POPULATION
		for (int j = 0; j < SIZE-2;j++)
		{
			childPopulation[i][j] = newChromosome[j];
		}
	}
	//3. MUTATE CHILD POPULATION
	//3a. FOR EVERY MEMBER IN CHILD POPULATION, DO THE FOLLOWING
	//3b. FOR EVERY NUMBER IN THE CHROMOSOME, SET PROBABILITY OF CHANGE TO 100
	//3c. IF NUMBER CHANGES, CHANGE TO NUMBER THAT MEETS DEGREE REQUIREMENTS
	//4. SET POPULATION TO CHILD POPULATION
	//5. GENERATION += 1, REPEAT WHILE LOOP
	
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

void displayArray(int arr[], int size)
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
