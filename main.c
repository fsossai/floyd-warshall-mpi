#include "pfw.h"
#include <time.h>

void free_all(float *d, int *pred, struct LinesCompData *lcd);

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&P);

	N = 0;
    if (IS_MASTER && argc == 2)
		N = atoi(argv[1]);

	MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	if (N == 0)
	{
		if (IS_MASTER)
		{
			fprintf(stderr, "MASTER: Wrong parameters\n");
			fprintf(stderr,"Usage: pfw <n>\t\tN: number of nodes in the graph\n");
		}
		
		MPI_Finalize();
		return 0;
	}

	float *d = NULL;			// vettore (matrice in row-major) delle distanze minime (allocato solo da MASTER)
	int *pred = NULL;			// vettore (matrice in row-major) del percorso minimo (allocato solo da MASTER),
								// pred(i,j) rappresenta il nodo precedente a 'j' nel percoso minimo da 'i' a 'j'
								// da questa matrice e' possibile ricavare i percorsi minimi passando da pred(i,j) a pred(i,pred(i,j)) ...
	int hprocess;				// ID del processo che possiede la versione aggiornata delle h-esime righe di 'd' e 'pred'
	int k;
	int hindex;					// riga di 'lcd.d' (e 'lcd.pred') che corrisponde alla h-esima riga in 'd' e 'pred'
	struct LinesCompData lcd;	// dati necessari al calcolo della partizione finale di 'd' e 'pred' da parte di uno SLAVE
	int* workload;				// workload[p] = numero di celle di 'd' e 'pred' che il processo 'p' deve calcolare
	int* displacement;			// displacement[p] = offset del primo byte della prima riga assegnata al proc. 'p' in 'd' e 'pred'
	
	workload = (int*)malloc(P * sizeof(int));
	displacement = (int*)malloc((P + 1) * sizeof(int));

	// Assegnamento dei carichi di lavoro
	for (int p = 0; p<P; p++)
		workload[p] = N * ( (N / P) + ( (p < N % P) ? 1 : 0 ) ); // workload[p] = N* ( partebassa(N/P) (+1 se p < N mod P) )
	displacement[0] = 0;
	for (int p = 1; p<=P; p++)
		displacement[p] = displacement[p-1] + workload[p-1];
	lcd.lines = workload[rank] / N;

	if (IS_MASTER)
	{
		printf("N = %i, P = %i : initializing ... \n", N, P);
		pred = (int*)malloc(N * N * sizeof(int));
		d = (float*)malloc(N * N * sizeof(float));
		random_mat(d,N,10);
		for (int i = 0; i<N; i++) // inizializzazione standard di 'pred': 'i' precede 'j' nel percorso i=>j
			for (int j = 0; j<N; j++)
				pred[INDEX(i,j)] = i;
	}

	init_lcd(&lcd,workload[rank]);
	MPI_Barrier(MPI_COMM_WORLD);

	if (IS_MASTER)
		printf("Computing ...\n");

	clock_t t = clock();
	double mpiTime = - MPI_Wtime();

	// MASTER puo' gia' salvarsi una copia della h-esima riga di 'd' e 'pred'
	if (IS_MASTER)		
	{
		memcpy(lcd.d_h, &d[INDEX(0,0)], N * sizeof(float));
		memcpy(lcd.pred_h, &pred[INDEX(0,0)], N * sizeof(float));
	}

	// MASTER distribuisce le partizioni delle tabelle iniziali
	MPI_Scatterv(d, workload, displacement,
					MPI_FLOAT, lcd.d, workload[rank],
					MPI_FLOAT, MASTER, MPI_COMM_WORLD);
	MPI_Scatterv(pred, workload, displacement,
				MPI_INT, lcd.pred, workload[rank],
				MPI_INT, MASTER, MPI_COMM_WORLD);

	for (int h = 0, k = 1; h<N; h++)
	{
		if (h * N < displacement[k]) 	// controllo quale processo possiede il più recente calcolo delle righe 'h' di 'd' e 'pred'
			hprocess = k-1;
		else
			hprocess = k++;
		
		// chi possiede la riga 'h' di 'd' e 'pred' la manda in broadcast agli altri processi
		if (rank == hprocess)
		{
			hindex = h - displacement[hprocess] / N;
			memcpy(lcd.d_h, &lcd.d[INDEX(hindex,0)], N * sizeof(float));
			memcpy(lcd.pred_h, &lcd.pred[INDEX(hindex,0)], N * sizeof(int));
		}
		MPI_Bcast(lcd.d_h, N, MPI_FLOAT, hprocess, MPI_COMM_WORLD);
		MPI_Bcast(lcd.pred_h, N, MPI_INT, hprocess, MPI_COMM_WORLD);

		compute_pfw(&lcd, h); 	// h-esima operazione triangolare
	}

	// MASTER raccoglie le partizioni finali
	MPI_Gatherv(lcd.d, workload[rank], MPI_FLOAT,
			d, workload, displacement, MPI_FLOAT,
			MASTER, MPI_COMM_WORLD);
	MPI_Gatherv(lcd.pred, workload[rank], MPI_INT,
			pred, workload, displacement, MPI_INT,
			MASTER, MPI_COMM_WORLD);

	mpiTime += MPI_Wtime();
	t = clock() - t;

	if (IS_MASTER)
		printf("P0 (MASTER):\n CPU time\t: %f ms\n MPI time\t: %f ms\n",
			1.e3 * ((double)t / CLOCKS_PER_SEC),
			mpiTime * 1.e3);

	free(workload);
	free(displacement);
	free_all(d, pred, &lcd);

	MPI_Finalize();
	return 0;
}

void free_all(float *d, int *pred, struct LinesCompData *lcd)
{
	if (IS_MASTER)
	{
		free(d);
		free(pred);
	}
	free(lcd->d);
	free(lcd->pred);
	free(lcd->d_h);
	free(lcd->pred_h);
}