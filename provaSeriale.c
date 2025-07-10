#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

// Genera un numero random tra min e max
double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// Genera un punto casuale nell'ipercubo k-dimensionale
Point* generate_random_point(Bounds *bounds, int k) {
    Point *p = create_point(k);
    for (int i = 0; i < k; i++) {
        p->coords[i] = random_double(bounds[i].min, bounds[i].max);
    }
    return p;
}


// Calcola il volume dell'ipercubo k-dimensionale
double calculate_hypercube_volume(Bounds *bounds, int k) {
    double volume = 1.0;
    for (int i = 0; i < k; i++) {
        volume *= (bounds[i].max - bounds[i].min);
    }
    return volume;
}

// ESEMPI DI FUNZIONI FORMA - Modifica queste per la tua figura specifica

// Esempio 1: Ipersfera k-dimensionale centrata nell'origine
int is_inside_hypersphere(Point *p, double radius) {
    double sum_squares = 0.0;
    for (int i = 0; i < p->dimension; i++) {
        sum_squares += p->coords[i] * p->coords[i];
    }
    return sum_squares <= radius * radius;
}

// Esempio 2: Ellissoide k-dimensionale
int is_inside_ellipsoid(Point *p, double *semi_axes) {
    double sum = 0.0;
    for (int i = 0; i < p->dimension; i++) {
        double term = p->coords[i] / semi_axes[i];
        sum += term * term;
    }
    return sum <= 1.0;
}

/* Esempio 3: Funzione generica per forme più complesse
// Puoi modificare questa per la tua figura specifica
int is_inside_custom_shape(Point *p) {
    // Esempio: una forma definita da una disequazione
    // Modifica questa funzione per la tua figura
    
    // Esempio semplice: un "diamante" k-dimensionale
    double sum_abs = 0.0;
    for (int i = 0; i < p->dimension; i++) {
        sum_abs += fabs(p->coords[i]);
    }
    return sum_abs <= 1.0;
}
    */

// Funzione principale per il metodo Monte Carlo
double monte_carlo_volume(Bounds *bounds, int k, long n_samples, 
                         int (*shape_function)(Point*)) {
    
    srand(time(NULL)); // Inizializza il generatore di numeri casuali
    
    long points_inside = 0;
    double hypercube_volume = calculate_hypercube_volume(bounds, k);
    
    printf("Inizio simulazione Monte Carlo...\n");
    printf("Dimensioni: %d\n", k);
    printf("Campioni: %ld\n", n_samples);
    printf("Volume ipercubo: %.6f\n\n", hypercube_volume);
    
    // Genera punti e conta quelli dentro la figura
    for (long i = 0; i < n_samples; i++) {
        Point *p = generate_random_point(bounds, k);
        
        if (shape_function(p)) {
            points_inside++;
        }
        
        free_point(p);
        
        // Progresso ogni 10% dei campioni
        if (i % (n_samples / 10) == 0 && i > 0) {
            double current_estimate = hypercube_volume * 
                                    ((double)points_inside / i);
            printf("Progresso: %.1f%% - Stima corrente: %.6f\n", 
                   (double)i / n_samples * 100, current_estimate);
        }
    }
    
    double estimated_volume = hypercube_volume * 
                             ((double)points_inside / n_samples);
    
    printf("\nRisultati finali:\n");
    printf("Punti totali: %ld\n", n_samples);
    printf("Punti dentro la figura: %ld\n", points_inside);
    printf("Rapporto: %.6f\n", (double)points_inside / n_samples);
    printf("Volume stimato: %.6f\n", estimated_volume);
    
    return estimated_volume;
}

// Funzione per testare la convergenza con diversi numeri di campioni
void test_convergence(Bounds *bounds, int k, 
                     int (*shape_function)(Point*)) {
    
    printf("\n=== TEST DI CONVERGENZA ===\n");
    long test_samples[] = {1000, 10000, 100000, 1000000};
    int n_tests = sizeof(test_samples) / sizeof(test_samples[0]);
    
    for (int i = 0; i < n_tests; i++) {
        printf("\n--- Test con %ld campioni ---\n", test_samples[i]);
        double volume = monte_carlo_volume(bounds, k, test_samples[i], 
                                         shape_function);
    }
}

// Esempio di uso per ipersfera
double hypersphere_wrapper(Point *p) {
    return is_inside_hypersphere(p, 1.0); // Raggio = 1
}

/* Esempio di uso per ellissoide
double *global_semi_axes; // Variabile globale per l'esempio
double ellipsoid_wrapper(Point *p) {
    return is_inside_ellipsoid(p, global_semi_axes);
}
*/
int main() {
    printf("=== METODO MONTE CARLO PER VOLUMI K-DIMENSIONALI ===\n\n");
    
    // ESEMPIO 1: Ipersfera 3D con raggio 1
    printf("ESEMPIO 1: Ipersfera 3D (raggio = 1)\n");
    printf("Volume teorico: %.6f\n", 4.0/3.0 * M_PI); // 4/3 * π
    
    int k = 3;
    Bounds bounds_sphere[] = {{-1.1, 1.1}, {-1.1, 1.1}, {-1.1, 1.1}};
    
    double volume_sphere = monte_carlo_volume(bounds_sphere, k, 1000000, 
                                            hypersphere_wrapper);
    
    printf("Errore relativo: %.2f%%\n", 
           fabs(volume_sphere - 4.0/3.0*M_PI) / (4.0/3.0*M_PI) * 100);
    
    
    /* ESEMPIO 2: Ellissoide 2D
    printf("\n\nESEMPIO 2: Ellissoide 2D (semi-assi: 2, 1)\n");
    printf("Volume teorico: %.6f\n", M_PI * 2 * 1); // π * a * b
    
    k = 2;
    double semi_axes[] = {2.0, 1.0};
    global_semi_axes = semi_axes;
    Bounds bounds_ellipse[] = {{-2.1, 2.1}, {-1.1, 1.1}};
    
    double volume_ellipse = monte_carlo_volume(bounds_ellipse, k, 500000, 
                                             ellipsoid_wrapper);
    
    printf("Errore relativo: %.2f%%\n", 
           fabs(volume_ellipse - M_PI*2*1) / (M_PI*2*1) * 100);
    
    
    // ESEMPIO 3: Forma personalizzata (diamante k-dimensionale)
    printf("\n\nESEMPIO 3: Diamante 4D\n");
    k = 4;
    Bounds bounds_diamond[] = {{-1.1, 1.1}, {-1.1, 1.1}, 
                              {-1.1, 1.1}, {-1.1, 1.1}};
    
    double volume_diamond = monte_carlo_volume(bounds_diamond, k, 500000, 
                                             is_inside_custom_shape);
    
    // Test di convergenza per la forma personalizzata
    test_convergence(bounds_diamond, k, is_inside_custom_shape);
    */
    return 0;
}