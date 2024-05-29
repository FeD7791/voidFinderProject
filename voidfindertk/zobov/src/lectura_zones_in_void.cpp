#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main() {
    // Nombre del archivo de entrada
    string inputFileName = "out_zones_in_void.dat";
    // Nombre del archivo de salida en formato ASCII
    string outputFileName = "txt_out_zones_in_void.txt";

    // Abrir el archivo de entrada en modo binario
    ifstream inputFile(inputFileName, ios::binary);
    if (!inputFile) {
        cerr << "Error al abrir el archivo de entrada: " << inputFileName << endl;
        return 1;
    }

    // Leer el número total de zonas
    int nzones;
    inputFile.read(reinterpret_cast<char*>(&nzones), sizeof(int));

    // Vector para almacenar las zonas y sus valores de nhl
    vector<pair<int, vector<int>>> zones(nzones);

    // Leer las zonas
    for (int h = 0; h < nzones; ++h) {
        // Leer el número total de elementos en la zona
        int nhl;
        inputFile.read(reinterpret_cast<char*>(&nhl), sizeof(int));

        // Vector para almacenar los elementos de la zona actual
        vector<int> zone(nhl);

        // Leer los elementos de la zona actual
        inputFile.read(reinterpret_cast<char*>(zone.data()), nhl * sizeof(int));

        // Guardar el par (nhl, zona) en el vector de zonas
        zones[h] = make_pair(nhl, zone);
    }

    // Cerrar el archivo de entrada
    inputFile.close();

    // Abrir el archivo de salida en modo texto
    ofstream outputFile(outputFileName);
    if (!outputFile) {
        cerr << "Error al abrir el archivo de salida: " << outputFileName << endl;
        return 1;
    }
    outputFile << "\n" << " ";
    outputFile << " nzones "<<nzones << " ";
    outputFile << "\n" << " ";

    // Escribir las zonas en el archivo de salida en formato ASCII
    int void_numero=0;
    for (const auto& zone : zones) {
        // Escribir el valor de nhl
        outputFile << " id void " << void_numero<< " zonas ";
        //outputFile << "voids  \n" << " ";

        //outputFile << zone.first << " ";
        //outputFile << "\n" << " ";
        //outputFile << "zona list  \n" << " ";

        // Escribir los elementos de la zona
        for (int element : zone.second) {
            outputFile << element << " " ;
        }
        outputFile << endl;
        void_numero = void_numero +1;
    }

    // Cerrar el archivo de salida
    outputFile.close();

    cout << "Datos guardados correctamente en " << outputFileName << endl;

    return 0;
}


