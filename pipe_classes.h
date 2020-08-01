#pragma once
#include<iostream>
using namespace std;

class Node
{
private:
	int num;
	double X, Y, Q;

public:
	Node()           // default constructor
	{
		num = 0;
		X = 0;
		Y = 0;
		Q = 0;
	}

	Node(int nodenum, double xvalue, double yvalue, double nodeload);
	double getX();    // gets X coordinate of the node created
	double getY();
	int getnodeno();    //gets the node number in the input file
	double getQ();	
};


class Tube
{
private:
	int tubenum;
	Node* nodeA;
	Node* nodeB;
	double d, q;

public:
	Tube()                    // default constructor
	{
		tubenum = 0;
		nodeA = NULL;
		nodeB = NULL;
	}
	Tube(int Tnum, Node* n1, Node* n2, double dia);
	double calcLength();
	double getd();
	double calcB();          // calculted the B value of the tube
	int getnode1();          // gets first node number of the tube
	int getnode2();
};


class PipeNet
{

private:
	Node * * vec_nodes;
	Tube** vec_tubes;
	int n_nodes;
	int n_tubes;

public:
	PipeNet(Node narray[], Tube tarray[], int nodes, int tubes);
	~PipeNet();          // destructor
	void calcFlux();      // calculates the q value in all the tubes
	
};