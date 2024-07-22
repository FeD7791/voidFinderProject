#include <iostream>
#include <fstream>
#include <vector>

struct Zone {
    int np;
    std::vector<int> m;
};

int main() {
    // Abre el archivo generado por el código C++
    std::ifstream input("out_particle_zone.dat", std::ios::binary);
    if (!input) {
        std::cerr << "No se pudo abrir el archivo de entrada.\n";
        return 1;
    }

    // Lee los datos desde el archivo binario
    int np, nzones;
    input.read(reinterpret_cast<char*>(&np), sizeof(int));
    input.read(reinterpret_cast<char*>(&nzones), sizeof(int));

    std::vector<Zone> zones(nzones);
    for (int h = 0; h < nzones; ++h) {
        input.read(reinterpret_cast<char*>(&zones[h].np), sizeof(int));
        zones[h].m.resize(zones[h].np);
        input.read(reinterpret_cast<char*>(zones[h].m.data()), sizeof(int) * zones[h].np);
    }

    input.close();

    // Abre un archivo de salida en formato ASCII
    std::ofstream output("txt_out_particle_zone2.txt");
    if (!output) {
        std::cerr << "No se pudo abrir el archivo de salida.\n";
        return 1;
    }

    // Escribe los datos en el archivo de salida en formato ASCII
    output << "\n np"<< std::endl; 
    output << np << std::endl;
    output << "\n nzones"<< std::endl; 
    output << nzones << std::endl;
    output << "\n"<< std::endl;

    for (const auto& zone : zones) {
        output << "\n------------------------"<< std::endl;
        output << " zone "<<zone.np << std::endl;
        output << " particulas  "<<std::endl;

        for (int value : zone.m) {
            output << value <<" ";
            //<< std::endl;
        }
    }

    output.close();

    std::cout << "Datos guardados en txt_out_particle_zone.txt\n";
    return 0;
}