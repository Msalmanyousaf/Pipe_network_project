#include<iostream>
#include<cmath>
#include"pipe_classes.h"
using namespace std;


Node::Node(int nodenum, double xvalue, double yvalue, double nodeload)
{
	num = nodenum;
	X = xvalue;
	Y = yvalue;
	Q = nodeload;
}

double Node:: getX()
{
	return X;
}

double Node::getY()
{
	return Y;
}

int Node::getnodeno()
{
	return num;
}

double Node::getQ()
{
	return Q;
}



Tube::Tube(int Tnum, Node* n1, Node* n2, double dia)
{
	tubenum = Tnum;
	nodeA = n1;
	nodeB = n2;
	d = dia;
}

double Tube::calcLength()
{
	double value1 = (nodeA->getX() - nodeB->getX());
	double value2 = (nodeA->getY() - nodeB->getY());
	return sqrt(value1*value1 + value2 * value2);
}

double Tube::getd()
{
	return d;
}

double Tube::calcB()
{
	double B;
	B = (3.14*9.81*pow(d, 4)) / (128 * 1E-6*calcLength());
	return B;
}

int Tube::getnode1()
{
	return nodeA->getnodeno();
}

int Tube::getnode2()
{
	return nodeB->getnodeno();
}



PipeNet::PipeNet(Node narray[], Tube tarray[], int nodes, int tubes)
{
	n_nodes = nodes;
	n_tubes = tubes;
	vec_nodes = new Node*[n_nodes];
	vec_tubes = new Tube*[n_tubes];
	for (int i = 0; i < n_nodes; i++)
	{
		vec_nodes[i] = &narray[i];
	}
	for (int i = 0; i < n_tubes; i++)
	{
		vec_tubes[i] = &tarray[i];
	}
}

PipeNet::~PipeNet()
{
	delete[] vec_nodes;
	delete[] vec_tubes;
}


void PipeNet::calcFlux()
{
	double* Bi = new double[n_tubes];
	for (int i = 0; i < n_tubes; i++)
	{
		Bi[i] = vec_tubes[i]->calcB();

	}
	double** B;
	B = new double*[n_nodes];          // B is node number by node number matrix
	for (int i = 0; i < n_nodes; i++)
	{
		B[i] = new double[n_nodes];
	}
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < n_nodes; j++)
		{
			B[i][j] = 0;
		}
	}
	for (int i = 0; i < n_tubes; i++)          // assembles B matrix
	{
		int n_a = vec_tubes[i]->getnode1();    
		int n_b = vec_tubes[i]->getnode2();
		B[n_a - 1][n_a - 1] += Bi[i];

		B[n_b - 1][n_b - 1] += Bi[i];
		B[n_a - 1][n_b - 1] -= Bi[i];
		B[n_b - 1][n_a - 1] -= Bi[i];
	}
	double* Qi = new double[n_nodes];        // Qi is load vector
	for (int i = 0; i < n_nodes; i++)
	{
		Qi[i] = vec_nodes[i]->getQ();
	}
	//boundary conditions applications.
	for (int i = 1; i < n_nodes; i++)
	{
		B[i][0] = 0;
		B[0][i] = 0;
	}
	B[0][0] = 1;
	Qi[0] = 0;
	double* minusQi = new double[n_nodes];
	for (int i = 0; i < n_nodes; i++)
	{
		minusQi[i] = -Qi[i];
	}
	// Gauss elimination on Bh+q=0 or Bh=minusQ

	// LU decomposition without pivoting 
	for (int k = 0; k < n_nodes - 1; k++)
	{

		for (int i = k + 1; i < n_nodes; i++)
		{
			if (B[i][k] != 0)
			{
				double mult = B[i][k] / B[k][k];
				B[i][k] = mult;
				for (int j = k + 1; j < n_nodes; j++)
					B[i][j] -= mult * B[k][j];
			}
		}
	}
	// forwad substitution for L y = b. y still stored in bb
	for (int i = 1; i < n_nodes; i++)
		for (int j = 0; j < i; j++) minusQi[i] -= B[i][j] * minusQi[j];
	// back substitution for U x = y. x still stored in bb
	for (int i = n_nodes - 1; i >= 0; i--)
	{
		for (int j = i + 1; j<n_nodes; j++) minusQi[i] -= B[i][j] * minusQi[j];
		minusQi[i] /= B[i][i];
	}

	cout << "\nNo of nodes are:\t" << n_nodes;
	cout << "\nNo of Tubes are:\t" << n_tubes;
	cout << "\nLoad vector (flow rates at nodes are)\n";
	for (int i = 0; i < n_nodes; i++)
	{
		cout << "At node" << i + 1 << "\t" << vec_nodes[i]->getQ() << " m3/s" << endl;
	}
	cout << "\nh values are\n";
	for (int i = 0; i < n_nodes; i++)
	{
		cout << "h" << i + 1 << "=\t" << minusQi[i] << " m" << endl;

	}
	double* q = new double[n_tubes];
	for (int i = 0; i < n_tubes; i++)
	{
		int nodeno1 = vec_tubes[i]->getnode1();
		int nodeno2 = vec_tubes[i]->getnode2();
		q[i] = Bi[i] * (minusQi[nodeno1 - 1] - minusQi[nodeno2 - 1]);
	}
	cout << "\nvalues of tube flows q\n";
	for (int i = 0; i < n_tubes; i++)
	{
		cout << "q" << i + 1 << "=\t" << q[i] << " m3/s" << endl;
	}
	delete[] q;
	delete[] minusQi;
	delete[] Qi;
	for (int i = 0; i < n_nodes; i++)
		delete[]B[i];

	delete[] B;
	delete[]Bi;
}
