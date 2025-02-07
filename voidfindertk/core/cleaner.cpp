#include "Catalogue.h"
#include <vector>
#include <string>

using namespace cbl;
using namespace catalogue;

extern "C" {

void process_catalogues(
    const char* file_voids,
    const char* file_tracers,
    double ratio,
    bool initial_radius,
    const double* delta_r_min,
    const double* delta_r_max,
    double threshold,
    const char* output_path,
    bool ol_crit,
    bool rescale,
    bool checkoverlap
    ) {

    try {
        std::vector<double> delta_r_vec = {*delta_r_min, *delta_r_max};
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

        if (ol_crit){
            void_catalogue.clean_void_catalogue(initial_radius, delta_r_vec, threshold, rescale, input_tracersCata, ChM, ratio, checkoverlap, cbl::catalogue::Var::_DensityContrast_);
            var_names_voids.emplace_back(cbl::catalogue::Var::_DensityContrast_);
        }else{
            void_catalogue.clean_void_catalogue(initial_radius, delta_r_vec, threshold, rescale, input_tracersCata, ChM, ratio, checkoverlap, cbl::catalogue::Var::_CentralDensity_);
            var_names_voids.emplace_back(cbl::catalogue::Var::_CentralDensity_);
        }
        //Print the received file paths

        // About clean_void_catalogue --cbl version 1-- Definition at line 1328 of file VoidCatalogue.cpp.
        // Parameters
        // initial_radius	erase voids outside a given interval delta_r of initial radius;
        // delta_r	the interval of accepted radii
        // threshold	the density threshold
        // rescale	true = for each void finds the larger radius enclosing density = threshold, false = skip the step
        // tracers_catalogue	object of class Catalogue with the tracers defining the void distribution (necessary if rescale = true)
        // ChM	object of ChainMesh3D class
        // ratio	distance from the void centre at which the density contrast is evaluated in units of the void radius. Ex: ratio = 0.1 \(\rightarrow\) 10% of the void radius lenght
        // checkoverlap	true \(\rightarrow\) erase all the voids wrt a given criterion, false \(\rightarrow\) skip the step
        // ol_criterion	the criterion for the overlap step (valid criteria: Var::DensityContrast, Var::CentralDensity)

        // Save the catalogue data to the provided output path
        void_catalogue.write_data(output_path, var_names_voids);

    }catch (cbl::glob::Exception &exc) {
        std::cerr << exc.what() << std::endl;
        exit(1);
    } catch (const std::exception &exc) {
    std::cerr << "Error: " << exc.what() << std::endl;
    exit(1);
    }
}

}

