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
    int np, nzones;
    fread(&np, sizeof(int), 1, input);
    fread(&nzones, sizeof(int), 1, input);

    struct Zone *zones = (struct Zone*)malloc(nzones * sizeof(struct Zone));
    for (int h = 0; h < nzones; ++h) {
        fread(&zones[h].np, sizeof(int), 1, input);
        zones[h].m = (int*)malloc(zones[h].np * sizeof(int));
        fread(zones[h].m, sizeof(int), zones[h].np, input);
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
    fprintf(output, "\n nzones\n");
    fprintf(output, "%d\n", nzones);
    fprintf(output, "\n");

    for (int i = 0; i < nzones; ++i) {
        fprintf(output, "\n------------------------\n");
        fprintf(output, " zone %d\n", zones[i].np);
        fprintf(output, " particulas  \n");

        for (int j = 0; j < zones[i].np; ++j) {
            fprintf(output, "%d ", zones[i].m[j]);
        }
        fprintf(output, "\n");
    }

    fclose(output);

    printf("Datos guardados en %s\n", outputFileName);

    // Liberar memoria
    for (int i = 0; i < nzones; ++i) {
        free(zones[i].m);
    }
    free(zones);

    return 0;
}

// int main() {
//     const char* inputFileName = "out_particle_zone.dat";
//     const char* outputFileName = "txt_out_particle_zone2.txt";
    
//     if (convertBinaryToAscii(inputFileName, outputFileName) != 0) {
//         return 1;
//     }

//     return 0;
// }