#include<iostream>
#include<fstream>
#include"Class_functions_definitions.h"

using namespace std;
int main()
{
	fstream mainfile;
	int n, t;                                // number of tubes and nodes
	double* mainarray;
	mainfile.open("pipedata.txt");         // read n and t from input file
	mainfile >> n;
	mainfile >> t;
	mainfile.close();

	fstream mainfile1;
	mainfile1.open("pipedata.txt");      // read pipe flow and pipe connectivity data from input file
	mainarray = new double[(n) * 3 + t * 3 + 2];            // main array to hold the data from file
	for (int i = 0; i < (n * 3 + t * 3 + 2); i++)
	{
		mainfile1 >> mainarray[i];
	}
	mainfile1.close();

	Node* nodearray = new Node[n];         // objects of Node class are created in the form of array
	Tube* tubearray = new Tube[t];          // objects of Tube class are created in the form of array

	int counter = 0;
	for (int i = 1; i < n + 1; i++)
	{
		double xcord = mainarray[counter + 2];
		double ycord = mainarray[counter + 2 + 1];
		double nodeflow = mainarray[counter + 2 + 1 + 1];

		Node a(i, xcord, ycord, nodeflow);
		nodearray[i - 1] = a;
		counter = counter + 3;
	}

	int mycounter = 0;
	for (int c = 0; c <t; c++)
	{
		int rand1 = mainarray[3 * n + mycounter + 2];             // Tube first node number
		int rand2 = mainarray[3 * n + mycounter + 1 + 2];         // Tube second node number
		double rand3 = mainarray[3 * n + mycounter + 2 + 2];       // Tube diameter
		Tube tt(c + 1, &nodearray[rand1 - 1], &nodearray[rand2 - 1], rand3);     // c+1 is tube number
		tubearray[c] = tt;
		mycounter = mycounter + 3;
	}

	PipeNet givennetwerk(nodearray, tubearray, n, t);
	givennetwerk.calcFlux();
	
	delete[] nodearray;
	delete[] tubearray;
	delete[] mainarray;
	system("pause");
	return 0;
}
