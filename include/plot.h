#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <set>
#include "purkinje.h"

using namespace std;

void printCellsVTK (Graph *g, int exec_number);
void printInfoModel (Graph *g, double *y0, double Dx, double Dt, double t_max, int exec_number, int num_fibers, double probGap, double d_bundle, double l_cell, double d_cell, double fl, double s_c, double s_p, double s_t, double C);
void writeGraphic (Graph *g, double t, Node **plot, int num_eq, int exec_number);
void writeIteration (Graph *g, int j, int M, int exec_number);
void writeVTK (Graph *g, int j, int exec_number);
void makePlot (Node **plot, int exec_number);
void printCellsRetroprogVTK (Graph *g, int exec_number);
void writeCylinderVTK (Graph *g, int exec_number);
void writeTime (clock_t begin, clock_t end, int exec_number);
void writeVelocity (double velocity, double delta_x, double delta_t, double t1, double t2, double *maxV, int exec_number);
void printExecution (int i, bool retro);
void printCooldown (Graph *g, int exec_number);

#endif
