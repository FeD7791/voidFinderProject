#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int first;
    int *second;
} pair;

void process_files(const char *inputFileName, const char *outputFileName) {
    // Abrir el archivo de entrada en modo binario
    FILE *inputFile = fopen(inputFileName, "rb");
    if (!inputFile) {
        fprintf(stderr, "Error al abrir el archivo de entrada: %s\n", inputFileName);
        return;
    }

    // Leer el número total de zonas
    int nzones;
    fread(&nzones, sizeof(int), 1, inputFile);

    // Vector para almacenar las zonas y sus valores de nhl
    pair zones[nzones];

    // Leer las zonas
    for (int h = 0; h < nzones; ++h) {
        // Leer el número total de elementos en la zona
        int nhl;
        fread(&nhl, sizeof(int), 1, inputFile);

        // Vector para almacenar los elementos de la zona actual
        int *zone = (int *)malloc(nhl * sizeof(int));

        // Leer los elementos de la zona actual
        fread(zone, sizeof(int), nhl, inputFile);

        // Guardar el par (nhl, zona) en el vector de zonas
        zones[h].first = nhl;
        zones[h].second = zone;
    }

    // Cerrar el archivo de entrada
    fclose(inputFile);

    // Abrir el archivo de salida en modo texto
    FILE *outputFile = fopen(outputFileName, "w");
    if (!outputFile) {
        fprintf(stderr, "Error al abrir el archivo de salida: %s\n", outputFileName);
        return;
    }
    fprintf(outputFile, "\n");
    fprintf(outputFile, " nzones %d ", nzones);
    fprintf(outputFile, "\n");

    // Escribir las zonas en el archivo de salida en formato ASCII
    int void_numero = 0;
    for (int i = 0; i < nzones; ++i) {
        // Escribir el valor de nhl
        //fprintf(outputFile, " id void %d zonas ", void_numero);
        fprintf(outputFile, "%d ", void_numero);

        // Escribir los elementos de la zona
        for (int j = 0; j < zones[i].first; ++j) {
            fprintf(outputFile, "%d ", zones[i].second[j]);
        }
        fprintf(outputFile, "\n");
        void_numero++;
    }

    // Cerrar el archivo de salida
    fclose(outputFile);

    printf("Datos guardados correctamente en %s\n", outputFileName);
}

// Example of how to call the function
// int main() {
//     process_files("out_zones_in_void.dat", "txt_out_zones_in_void.txt");
//     return 0;
// }