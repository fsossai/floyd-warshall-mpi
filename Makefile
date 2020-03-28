# Set all variables appropriately 

CC = mpicc
INCLUDE = .
LIBS = .

pfw: main.c pfw.c
	$(CC) -O3 -I$(INCLUDE) -L$(LIBS) $^ -o pfw