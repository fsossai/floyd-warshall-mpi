#include "pfw.h"

// Floyd-Warshall sequenziale
void compute_fw(float *d, int *pred)
{	
	for (int h = 0; h < N; h++)
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if ( d[INDEX(i,h)] + d[INDEX(h,j)] < d[INDEX(i,j)] )
				{
					d[INDEX(i,j)] = d[INDEX(i,h)] + d[INDEX(h,j)];
					pred[INDEX(i,j)] = pred[INDEX(h,j)];
				}
}

// Floyd-Warshall per il calcolo della h-esima operazione triangolare
void compute_pfw(struct LinesCompData *lcd, int h)
{
	for (int i = 0; i < lcd->lines; i++)
		for (int j = 0; j < N; j++)
			if ( lcd->d[INDEX(i,h)] + lcd->d_h[j] < lcd->d[INDEX(i,j)] )
			{
				lcd->d[INDEX(i,j)] = lcd->d[INDEX(i,h)] + lcd->d_h[j];
				lcd->pred[INDEX(i,j)] = lcd->pred_h[j];
			}
}

void init_lcd(struct LinesCompData *lcd, int workload)
{
	lcd->d = (float*)malloc(workload * sizeof(float));
	lcd->pred = (int*)malloc(workload * sizeof(int));
	lcd->d_h = (float*)malloc(N * sizeof(float));
	lcd->pred_h = (int*)malloc(N * sizeof(int));
}

// Generazione di un problema casuale con pesi sugli archi tra 0 e 'max' escluso
void random_mat(float *c, int N, int max)
{
	srand( time(0) );
	for (int i = 0; i<N; i++)
		for (int j = 0; j<N; j++)
			c[INDEX(i,j)] = rand() % max;
}

// Per la visualizzazione di una matrice di float
void show_mat_float(float *mat)
{
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
			printf("%.0f\t", mat[INDEX(i,j)]);
		printf("\n");
	}
}

// Per la visualizzazione di una matrice intera
void show_mat_int(int *mat)
{
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
			printf("%i\t", mat[INDEX(i,j)]);
		printf("\n");
	}
}