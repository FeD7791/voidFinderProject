#include "Catalogue.h"
#include <vector>
#include <string>

using namespace cbl;
using namespace catalogue;

extern "C" {

void process_catalogues(const char* file_voids, const char* file_tracers, 
                        double ratio, bool initial_radius, const double* delta_r, int delta_r_size, 
                        double threshold, const char* output_path) {

    try {
        std::vector<double> delta_r_vec(delta_r, delta_r + delta_r_size);

        // Load the input void catalogue
        std::vector<cbl::catalogue::Var> var_names_voids = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_Radius_};
        std::vector<int> columns_voids = {1, 2, 3, 4};
        cbl::catalogue::Catalogue void_catalogue = Catalogue(cbl::catalogue::ObjectType::_Void_, cbl::CoordinateType::_comoving_, var_names_voids, columns_voids, {file_voids}, 0);

        // Build the tracer catalogue
        std::vector<cbl::catalogue::Var> var_names_tracers = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_};
        std::vector<int> columns_tracers = {1, 2, 3};
        cbl::catalogue::Catalogue tracers_catalogue = Catalogue(cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_tracers}, 0);

        double mps = tracers_catalogue.mps();
        cbl::chainmesh::ChainMesh3D ChM(2*mps, tracers_catalogue.var(cbl::catalogue::Var::_X_), tracers_catalogue.var(cbl::catalogue::Var::_Y_), tracers_catalogue.var(cbl::catalogue::Var::_Z_), void_catalogue.Max(cbl::catalogue::Var::_Radius_));
        auto input_tracersCata = std::make_shared<cbl::catalogue::Catalogue>(cbl::catalogue::Catalogue(std::move(tracers_catalogue)));

        void_catalogue.clean_void_catalogue(initial_radius, delta_r_vec, threshold, true, input_tracersCata, ChM, ratio, true, cbl::catalogue::Var::_CentralDensity_);

        var_names_voids.emplace_back(cbl::catalogue::Var::_CentralDensity_);
        
        // Save the catalogue data to the provided output path
        void_catalogue.write_data(output_path, var_names_voids);

    } catch (cbl::glob::Exception &exc) {
        std::cerr << exc.what() << std::endl;
        exit(1);
    }
}

}





















// #include "Catalogue.h"
// #include <vector>
// #include <string>

// using namespace cbl;
// using namespace catalogue;

// extern "C" {

// void process_catalogues(const char* file_voids, const char* file_tracers, 
//                         double ratio, bool initial_radius, const double* delta_r, int delta_r_size, double threshold) {

//     try {
//         std::vector<double> delta_r_vec(delta_r, delta_r + delta_r_size);

//         // Load the input void catalogue
//         std::vector<cbl::catalogue::Var> var_names_voids = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_Radius_};
//         std::vector<int> columns_voids = {1, 2, 3, 4};
//         cbl::catalogue::Catalogue void_catalogue = Catalogue(cbl::catalogue::ObjectType::_Void_, cbl::CoordinateType::_comoving_, var_names_voids, columns_voids, {file_voids}, 0);

//         // Build the tracer catalogue
//         std::vector<cbl::catalogue::Var> var_names_tracers = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_};
//         std::vector<int> columns_tracers = {1, 2, 3};
//         cbl::catalogue::Catalogue tracers_catalogue = Catalogue(cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_tracers}, 0);

//         double mps = tracers_catalogue.mps();
//         cbl::chainmesh::ChainMesh3D ChM(2*mps, tracers_catalogue.var(cbl::catalogue::Var::_X_), tracers_catalogue.var(cbl::catalogue::Var::_Y_), tracers_catalogue.var(cbl::catalogue::Var::_Z_), void_catalogue.Max(cbl::catalogue::Var::_Radius_));
//         auto input_tracersCata = std::make_shared<cbl::catalogue::Catalogue>(cbl::catalogue::Catalogue(std::move(tracers_catalogue)));

//         void_catalogue.clean_void_catalogue(initial_radius, delta_r_vec, threshold, true, input_tracersCata, ChM, ratio, true, cbl::catalogue::Var::_CentralDensity_);

//         var_names_voids.emplace_back(cbl::catalogue::Var::_CentralDensity_);
//         std::string mkdir = "mkdir -p ../output/";
//         if (system(mkdir.c_str())) {}

//         std::string cata_out = "../output/cleaned_void_catalogue.out";
//         void_catalogue.write_data(cata_out, var_names_voids);

//     } catch (cbl::glob::Exception &exc) {
//         std::cerr << exc.what() << std::endl;
//         exit(1);
//     }
// }

// }

























// #include "Catalogue.h"

// using namespace cbl;
// using namespace catalogue;

// extern "C"{

// void process_catalogues(const std::string& file_voids, const std::string& file_tracers,
//                         double ratio, bool initial_radius, const std::vector<double>& delta_r, double threshold) {

//   try {

//     // ------------------------------------------
//     // ----- load the input void catalogue ------
//     // ------------------------------------------

//     // std::vector containing the variable name list to read from file
//     std::vector<cbl::catalogue::Var> var_names_voids = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_Radius_};
    
//     // std::vector containing the columns corresponding to each attribute
//     std::vector<int> columns_voids = {1, 2, 3, 4};
    
//     // catalogue constructor
//     cbl::catalogue::Catalogue void_catalogue = Catalogue(cbl::catalogue::ObjectType::_Void_, cbl::CoordinateType::_comoving_, var_names_voids, columns_voids, {file_voids}, 0);
     
//     // ---------------------------------------
//     // ----- build the tracer catalogue ------
//     // ---------------------------------------
      
//     // std::vector containing the variable name list to read from file
//     std::vector<cbl::catalogue::Var> var_names_tracers = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_};
      
//     // std::vector containing the column corresponding to each attribute
//     std::vector<int> columns_tracers = {1, 2, 3};

//     // catalogue constructor
//     cbl::catalogue::Catalogue tracers_catalogue = Catalogue(cbl::catalogue::ObjectType::_Halo_, cbl::CoordinateType::_comoving_, var_names_tracers, columns_tracers, {file_tracers}, 0);

//     // store the mean particle separation of the simulation
//     double mps = tracers_catalogue.mps();      

//     // generate the chain mesh of the input tracer catalogue
//     cbl::chainmesh::ChainMesh3D ChM(2*mps, tracers_catalogue.var(cbl::catalogue::Var::_X_), tracers_catalogue.var(cbl::catalogue::Var::_Y_), tracers_catalogue.var(cbl::catalogue::Var::_Z_), void_catalogue.Max(cbl::catalogue::Var::_Radius_));

//     // make a shared pointer to tracers_catalogue
//     auto input_tracersCata = std::make_shared<cbl::catalogue::Catalogue>(cbl::catalogue::Catalogue(std::move(tracers_catalogue)));

//     // --------------------------------------------
//     // ----- build the cleaned void catalogue -----
//     // --------------------------------------------

//     // catalogue constructor
//     void_catalogue.clean_void_catalogue(initial_radius, delta_r, threshold, true, input_tracersCata, ChM, ratio, true, cbl::catalogue::Var::_CentralDensity_);

//     // store the obtained catalogue in an ASCII file
//     var_names_voids.emplace_back(cbl::catalogue::Var::_CentralDensity_);
//     std::string mkdir = "mkdir -p ../output/";
//     if (system(mkdir.c_str())) {}

//     std::string cata_out = "../output/cleaned_void_catalogue.out";
//     void_catalogue.write_data(cata_out, var_names_voids);
    
//   }
//   catch (cbl::glob::Exception &exc) { 
//     std::cerr << exc.what() << std::endl; 
//     exit(1); 
//   }
// }
// }

// int main() {
//   // Define the parameters
//   std::string file_voids = "../input/void_catalogue.txt";
//   std::string file_tracers = "../input/halo_catalogue.txt";
//   double ratio = 1.5;
//   bool initial_radius = true;
//   std::vector<double> delta_r = {17., 150.};
//   double threshold = 0.3;

//   // Call the function with the parameters
//   process_catalogues(file_voids, file_tracers, ratio, initial_radius, delta_r, threshold);

//   return 0;
// }