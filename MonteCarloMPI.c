#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


    /*  -------------------------------------------------
    -----------------------STRUCT------------------------
    ----------------------------------------------------- */

    // Struttura per definire i limiti di ogni dimensione
    typedef struct {
        double min;
        double max;
    } Bounds;

    // Struttura per un punto k-dimensionale
    typedef struct {
        double *coords;
        int dimension;
    } Point;

    //Struttura per trasmettere i parametri dell'ellissoide k-dimensionale
    typedef struct {
        int n;             // numero di semiassi
        double *semi_axes;    // array dei semiassi
    } EllipsoidParams;


    /*  -------------------------------------------------
    -----------------Metodi di utilità-------------------
    ----------------------------------------------------- */

    // Funzione per creare un punto k-dimensionale
    Point* create_point(int k) {
        Point *p = malloc(sizeof(Point));
        p->coords = malloc(k * sizeof(double));
        p->dimension = k;
        return p;
    }

    // Funzione per liberare la memoria di un punto
    void free_point(Point *p) {
        free(p->coords);
        free(p);
    }

        // Genera un numero random tra min e max (thread-safe)
    double random_double(double min, double max, unsigned int* seed) {
        double random_val = min + (max - min) * ((double)rand_r(seed) / RAND_MAX);
        // printf("Thread %d - Numero generato: %.6f\n", omp_get_thread_num(), random_val);
        return random_val;
    }
    // Genera un punto casuale nell'ipercubo k-dimensionale
    Point* generate_random_point(Bounds *bounds, int k, Point* p, unsigned int* seed) {
        for (int i = 0; i < k; i++) {
            p->coords[i] = random_double(bounds[i].min, bounds[i].max, seed);
        }
        return p;
    }

    Bounds* boundsBuilder(Bounds *bounds, int dimensione , double min , double max){
        for ( int i = 0 ; i< dimensione ; i++ ){
            bounds[i].min = min;
            bounds[i].max = max;
        }
        return bounds;
    }

    // Calcola il volume dell'ipercubo k-dimensionale
    double calculate_hypercube_volume(Bounds *bounds, int k) {
        double volume = 1.0;
        for (int i = 0; i < k; i++) {
            volume *= (bounds[i].max - bounds[i].min);
        }
        return volume;
    }

    //Ipersfera k-dimensionale 
    int is_inside_hypersphere(Point *p, double radius) {

        double sum_squares = 0.0;
        for (int i = 0; i < p->dimension; i++) {
            sum_squares += p->coords[i] * p->coords[i];
        }
        return sum_squares <= radius * radius;
    }

    //Ellissoide k-dimensionale
    int is_inside_ellipsoid(Point *p, int *semi_axes) {
        double sum = 0.0;
        for (int i = 0; i < p->dimension; i++) {
            double term = p->coords[i] / semi_axes[i];
            sum += term * term;
        }
        return sum <= 1.0;
    }

        // Esempio di uso per ipersfera
    double hypersphere_wrapper(Point *p, void* params) {
        double raggio = *((double*) params);
        return is_inside_hypersphere(p, raggio); // Raggio = 1
    }

    double ellipsoid_wrapper(Point *p, void* params) {
        EllipsoidParams* ellissoide = ((EllipsoidParams*) params);
        int* lunghezze_semiassi = (int*) malloc(sizeof(int)* ellissoide->n);
        for ( int i = 0 ;i < ellissoide->n ; i++  ){
            lunghezze_semiassi[i] = ellissoide->semi_axes[i];
        }
        return is_inside_ellipsoid(p, lunghezze_semiassi);
    }

    double computeSphereVolume(int dim , double raggio){
        return (pow(M_PI , dim/2.0) / tgamma(dim/2.0 + 1.0)) * pow(raggio , dim);

    }

    double computeEllipsoidVolume(int dim , double* array_of_semi_axes){
        double prodotto_semiassi = 1;    
        for ( int i = 0 ; i < dim ; i++ ){
                prodotto_semiassi *= array_of_semi_axes[i];
            }
        return (pow(M_PI, dim / 2.0) / tgamma(dim / 2.0 + 1.0)) * prodotto_semiassi;
    }

    double get_asse_max(double* array, int dimensione){

        double max = 0 ;
        for ( int i = 0 ; i < dimensione ; i++ ){
            if( array[i] > max ){ 
                max = array[i];
            }
        }
        return max;
    }

    void printStats( double volume_stimato , double volume_effettivo){
        printf("Errore relativo: %.2f%%\n", 
        fabs(volume_stimato - volume_effettivo) / (volume_effettivo) * 100);
    }



// Funzione che stima il volume di una figura in input
double monte_carlo_volume(Bounds *bounds, int k, long n_samples, 
                        double (*shape_function)(Point*, void*), void* params) {
    

    int rank,size;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
                            
    long sum_of_points_inside = 0;
    const double hypercube_volume = calculate_hypercube_volume(bounds, k);
    
    if ( 0 == rank){
        printf("Inizio simulazione Monte Carlo...con %ld campioni\n", n_samples);
        printf("Volume ipercubo: %.6f\n\n", hypercube_volume);
    }
        

    // Seed unico per ogni processo
    unsigned int seed = (unsigned int)(time(NULL) + rank * 997 + getpid());
    
    Point* p = create_point(k);
    long my_points_inside = 0;
    
    
    // Ogni processo esegue la sua porzione ( da modificare)
    long samples_per_thread = n_samples / size;
    long start_sample = rank * samples_per_thread;
    long end_sample = (rank == size - 1) ? n_samples : start_sample + samples_per_thread;
    
    for (long i = start_sample; i < end_sample; i++) {
        p = generate_random_point(bounds, k, p, &seed);
        if (shape_function(p, params)) {
            my_points_inside++;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&my_points_inside,       // Dato locale da inviare
                &sum_of_points_inside,  // Buffer che riceve il risultato (solo nel processo root)
                1,                      // Numero di elementi           
                MPI_LONG,                // Tipo di dato
                MPI_SUM,               //Tipo di operazione
                0,                      //Rank del processo che riceve il risultato ( root )
                MPI_COMM_WORLD          //Communicator (processi coinvolti)
            );
    
    free_point(p);

    
    

    const double estimated_volume = hypercube_volume * 
                            ((double)sum_of_points_inside / n_samples);
    if( rank == 0 ){
    printf("\nRisultati finali:\n");
    printf("Punti totali: %ld\n", n_samples);
    printf("Punti dentro la figura: %ld\n", sum_of_points_inside);
    printf("Volume stimato: %.6f\n", estimated_volume);
    
    }
    return estimated_volume;
}



    /*  -------------------------------------------------
    -----------------------MAIN--------------------------
    ----------------------------------------------------- */


    int main(int argc , char** argv) {
        
    int rank_main;
        MPI_Init(&argc, &argv);

        MPI_Comm_rank(MPI_COMM_WORLD, &rank_main);

        if ( rank_main == 0 ){
        printf("=== METODO MONTE CARLO PER VOLUMI N-DIMENSIONALI ===\n\n");
        /*Si passa da linea di comando 
            -il parametro dimensione (le dimensioni dello spazio vettoriale in cui lavoriamo)
            -il parametro samples, ovvero il numero di punti che vengono generati per calcolare il volume
            -il nome della figura di cui si vuole calcolare
            - i parametri della figura ( raggio  o lunghezza semi assi ) */
        printf("Variabili da linea di comando: samples , dim, figura , params\n");
        }
        int dim; // dimensione dello spazio vettoriale nel quale lavoriamo
        char figura[20];
        int samples;
        double tstart, tstop;
        if( argv[1] != NULL ){
             samples= atoi(argv[1]);
        }else{
            samples = 100000;
            if ( 0 == rank_main){
            printf("Samples non specificati, esecuzione con %d\n", samples);
            }    
        }
        if( argc > 2 &&  argv[2] != NULL ){
             dim= atoi(argv[2]);
        }else{
            dim = 3;
            if ( 0 == rank_main){    
            printf("Dimensione non specificata, esecuzione con dim = %d\n", dim);
            }
        }

        if( argc>3 && argv[3] != NULL){
            // Salvo la stringa passata in una variabile
            strncpy(figura, argv[3], sizeof(figura) - 1);
            figura[sizeof(figura) - 1] = '\0';  // assicuro la terminazione
        }else{
            strcpy(figura,"sfera");
            if ( 0 == rank_main){
            printf("Figura non specificata, esecuzione su %s\n", figura);
            }
        }
        
        Bounds *bounds =  (Bounds *) malloc(sizeof(Bounds)*dim);
        double volume_stimato;
        double volume_teorico;

        if( strcmp(figura, "sfera") == 0 && argv[4] != NULL ){
            
            //Acquisizione raggio, o inizializzazione "statica"
            double raggio;
            if ( argc < 4){
                raggio = 1;
                if ( 0 == rank_main){
                printf("Raggio non specificato, esecuzione con raggio = %f\n", raggio);
                }
            }else{
                raggio = atof(argv[4]);
            }
            
            //inizializzazione iper-cubo
            bounds = boundsBuilder( bounds, dim ,(-raggio)-1 , raggio +1 );
            if ( 0 == rank_main){
                printf("Ipersfera %dD (raggio = %f)\n", dim, raggio);
            }
            //calcolo volume teorico dell' ipersfera
            volume_teorico = computeSphereVolume(dim , raggio);
            if ( 0 == rank_main){
                printf("Volume teorico: %.6f\n", volume_teorico); 
            }
            //calcolo volume stimato con metodo di montecarlo
            tstart =MPI_Wtime();
            volume_stimato = monte_carlo_volume(bounds, dim, samples,hypersphere_wrapper, &raggio);
            tstop =MPI_Wtime();

            //stampa statistiche
            printStats(volume_stimato , volume_teorico);            
        }
        else if(strcmp(figura, "ellisse") == 0 && argv[4] != NULL ){
            //definizione parametri
            double asse_max;
            int numero_parametri_ellisse = argc - 4; // 4 sono i primi parametri: nome del programma, numero samples , dimensioni, numero della figura
            asse_max = 0 ;
            double * array_of_semi_axes = (double*) malloc(sizeof(double)*dim);
            
            //gestione dati incompleti
            if ( numero_parametri_ellisse != dim ){
                printf("Lunghezze dei semiassi non definite interamente o correttamente\n");
                free(array_of_semi_axes);
                return 1;
            }

            //acquisizione parametri semiassi
            for ( int i= 0 ; i < dim ; i++ ){
                 
                if ( atof(argv[4+i]) > 0 ){
                array_of_semi_axes[i] = atof(argv[4+i]);
                }
                else{
                    printf("I valori dei semiassi devono essere positivi!!\n");
                    free(array_of_semi_axes);
                    return 1;
                }
            }
            asse_max = get_asse_max(array_of_semi_axes, dim);

            //costruzione ipercubo
            bounds = boundsBuilder( bounds, dim, (- asse_max)-1 , asse_max +1);

            //inizializzazione struttura per i parametri dell'ellisse
            EllipsoidParams* ellissoide = (EllipsoidParams*) malloc(sizeof(EllipsoidParams));
            ellissoide->n = dim;
            ellissoide->semi_axes = array_of_semi_axes;
            if ( 0 == rank_main){
                printf("Elissoide %dD \n", ellissoide->n);
            }
            //calcolo volume ellissoide
            volume_teorico = computeEllipsoidVolume(dim , array_of_semi_axes);
            if ( 0 == rank_main){
                printf("Volume teorico: %.6f\n", volume_teorico); 
            }
            //calcolo volume stimato
            tstart = MPI_Wtime();
            volume_stimato = monte_carlo_volume(bounds, dim, samples, ellipsoid_wrapper, ellissoide);
            tstop =MPI_Wtime();
            

            //stampa statistiche
            if( 0 == rank_main ){
                printStats(volume_stimato , volume_teorico);
            }
            //libero la memoria
            free(array_of_semi_axes);
            free(ellissoide);
        }
        else{
            printf("%s non e' una figura ammessa", figura);
        }

        free(bounds);
        if ( 0 == rank_main){
            printf("Tempo di esecuzione: %.6f\n", tstop-tstart);
        }

        MPI_Finalize();
        return 0;
    }
