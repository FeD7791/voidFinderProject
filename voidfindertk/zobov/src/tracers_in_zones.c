#include <stdio.h>
#include <stdlib.h>



struct Zone {
    int np;
    int *m;
};

int get_tracers_in_zones(const char* inputFileName, const char* outputFileName) {
    // Abre el archivo generado por el código C++
    FILE *input = fopen(inputFileName, "rb");
    if (!input) {
        fprintf(stderr, "No se pudo abrir el archivo de entrada.\n");
        return 1;
    }

    // Lee los datos desde el archivo binario
    int np, n_npart;
    fread(&np, sizeof(int), 1, input);
    fread(&n_npart, sizeof(int), 1, input);

    struct Zone *npart = (struct Zone*)malloc(n_npart * sizeof(struct Zone));
    for (int h = 0; h < n_npart; ++h) {
        fread(&npart[h].np, sizeof(int), 1, input);
        npart[h].m = (int*)malloc(npart[h].np * sizeof(int));
        fread(npart[h].m, sizeof(int), npart[h].np, input);
    }

    fclose(input);

    // Abre un archivo de salida en formato ASCII
    FILE *output = fopen(outputFileName, "w");
    if (!output) {
        fprintf(stderr, "No se pudo abrir el archivo de salida.\n");
        return 1;
    }

    // Escribe los datos en el archivo de salida en formato ASCII
    fprintf(output, "\n np\n");
    fprintf(output, "%d\n", np);
    fprintf(output, "\n n_npart\n");
    fprintf(output, "%d\n", n_npart);
    fprintf(output, "\n");

    for (int i = 0; i < n_npart; ++i) {
        fprintf(output, "\n------------------------\n");
        fprintf(output, " Nparticles %d\n", npart[i].np);
        fprintf(output, " particulas  \n");

        for (int j = 0; j < npart[i].np; ++j) {
            fprintf(output, "%d ", npart[i].m[j]);
        }
        fprintf(output, "\n");
    }

    fclose(output);

    printf("Datos guardados en %s\n", outputFileName);

    // Liberar memoria
    for (int i = 0; i < n_npart; ++i) {
        free(npart[i].m);
    }
    free(npart);

    return 0;
}

// int main() {
//     const char* inputFileName = "part_vs_zone.dat";
//     const char* outputFileName = "txt_out_particle_zone2.txt";
    
//     if (get_tracers_in_zones(inputFileName, outputFileName) != 0) {
//         return 1;
//     }

//     return 0;
// }