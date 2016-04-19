#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <set>
#include "purkinje.h"

using namespace std;

void printCellsVTK (Graph *g, int exec_number);
void printInfoModel (Graph *g, double *y0, double Dx, double Dt, double t_max, int exec_number, int num_fibers, double probGap,
			double diameter_bundle, double s_c, double s_p, double s_t, int pac, int ret, int cool);
void writeGraphic (Graph *g, double t, int *id, int size, int num_eq, int exec_number);
void writeVTK (Graph *g, int j, int exec_number);
void makePlot (int *plot, int size, int exec_number);
void printCellsRetroprogVTK (Graph *g, int exec_number);
void writeCylinderVTK (Graph *g, int exec_number);
void writeVelocity (double velocity, double delta_x, double delta_t, double t1, double t2, double *maxV, int exec_number);
void printExecution (int i, bool retro);
void printCooldown (Graph *g, int exec_number);

#endif
