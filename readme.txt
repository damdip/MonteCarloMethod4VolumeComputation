
Monte Carlo Method for N-Dimensional Volume Computation
=======================================================

Questo progetto implementa un algoritmo generico per il calcolo del volume di figure n-dimensionali (es. ipersfera, ellissoide) usando il metodo Monte Carlo. Sono fornite due versioni parallele:

- MonteCarloOMP.c — versione con OpenMP (multi-threading su memoria condivisa)
- MonteCarloMPI.c — versione con MPI (calcolo distribuito su più processi)

-------------------------------------------------------
Compilazione
-------------------------------------------------------

OpenMP:
    gcc -fopenmp -O2 MonteCarloOMP.c -o mc_omp -lm

MPI:
    mpicc -O2 MonteCarloMPI.c -o mc_mpi -lm

-------------------------------------------------------
Esecuzione
-------------------------------------------------------

Esempio con OpenMP:
    OMP_NUM_THREADS=4 ./mc_omp 1000000 3 sfera

Esempio con MPI:
    mpirun -n 4 ./mc_mpi 1000000 3 sfera

Argomenti da riga di comando:

1. samples – numero di punti Monte Carlo da generare
2. dim     – dimensione dello spazio (es. 3 per 3D)
3. figura  – tipo di figura: "sfera" oppure "ellisse"
4. params  – eventuali parametri (es. semiassi per ellisse, opzionali)

-------------------------------------------------------
Funzionalità principali
-------------------------------------------------------

- Calcolo del volume tramite metodo Monte Carlo in R^n
- Supporto per figura "sfera" e "ellisse"
- Generalizzazione a dimensione arbitraria
- Supporto al parallelismo via OpenMP o MPI
- Calcolo del tempo di esecuzione
- Stampa di errore relativo rispetto al volume teorico (dove disponibile)

-------------------------------------------------------
Requisiti
-------------------------------------------------------

- GCC (con supporto a -fopenmp)
- MPI (es. OpenMPI, MPICH)
- Libreria math.h (link con -lm)

-------------------------------------------------------
File inclusi
-------------------------------------------------------

- MonteCarloOMP.c – versione OpenMP
- MonteCarloMPI.c – versione MPI
- readme.txt – questo file

-------------------------------------------------------
Autore
-------------------------------------------------------

Progetto didattico per lo studio del calcolo Monte Carlo e del parallelismo con MPI/OpenMP.


Generato con ChatGpt :)
