#ifndef MONODOMAIN_H
#define MONODOMAIN_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "purkinje.h"

using namespace std;

#define CYCLE_L 4000					// Intervalo de tempo do pacing
#define PI 3.14159265359			// Valor de pi

// Define o tipo Func, ponteiro para função que recebe t, e um vetor y como parâmetro e retorna a avaliação da função
typedef
	double (*Func) (int k, double t, double *y);

void solveEDO (Graph *g, double t, int k, Func *func, int num_eq);
void solveEDP (Graph *g, Func *func, double t, int k, int num_eq);
void makeMatrix_A (Graph *g);
void makeVector_b (Graph *g);
void solveEDP_Imp (Graph *g);
double* LU (int n);
void nextTimestep (Graph *g, int num_eq);
void setParameters (double Dx, double Dt, double cm);
void checkRetropropagation (Node *ptr);
bool printRetropropagation (Graph *g, int exec_number);
void printCooldown (Graph *g);

#endif
