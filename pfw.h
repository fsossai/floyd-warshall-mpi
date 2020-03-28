#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#define IS_MASTER !rank
#define IS_SLAVE rank
#define MASTER 0
#define INDEX(i,j) (i*N + j)

struct LinesCompData
{
	float *d_h;		// Memorizza l'h-esima riga della tabella delle distanze
	int *pred_h;	// Memorizza l'h-esima riga della tabella dei percorsi
	float *d;		// Vettore (matrice in row-major) della partizione della tablella delle distanze
	int *pred;		// Vettore (matrice in row-major) della partizione della tabella dei percorsi
	int lines;		// Numero di righe (di N elementi) contenute in 'd' e 'pred'
};

int P, rank, N;

void compute_fw(float *d, int *pred);
void compute_pfw(struct LinesCompData *lcd, int h);
void init_lcd(struct LinesCompData *lcd, int workload);
void random_mat(float *c, int N, int max);
void show_mat_float(float *mat);
void show_mat_int(int *mat);
