//********************************** VERSION 1.0 **********************************//
  
 void cbl::catalogue::Catalogue::clean_void_catalogue (const bool initial_radius, const std::vector<double> delta_r, const double threshold, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
 {
   const double start_time = omp_get_wtime();
     
   // ---------------------------------------------------- //
   // ---------------- Cleaning Procedure ---------------- //
   // ---------------------------------------------------- //
  
   cout << endl;
   coutCBL << par::col_blue << "> > > > > >>>>>>>>>> CLEANING PROCEDURE STARTED  <<<<<<<<< < < < < <" << par::col_default << endl << endl;
   coutCBL << "Voids in the initial Catalogue: " << nObjects() << endl << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
   coutCBL << par::col_green << " --- Removing spurious voids --- " << par::col_default << endl << endl;
   coutCBL << "Removed voids: " << endl;
  
   double spurious_time = omp_get_wtime();
   if (initial_radius) {
     vector<bool> remove(nObjects(), false);
     for (size_t i=0; i<nObjects(); i++) 
       if (radius(i) < delta_r[0] || delta_r[1] < radius(i)) remove[i] = true;
  
     remove_objects(remove);
     cout << "\t r_min-r_max criterion: " << count(remove.begin(), remove.end(), true) << endl;
   }
   
   compute_centralDensity(tracers_catalogue, ChM, threshold);
   coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
   coutCBL << "Time spent by the spurious voids-checking procedure: " << fixed << setprecision(3) << omp_get_wtime()-spurious_time << " seconds" << endl;
  
  
   // ---------------------------------------------------- //
   // ----------------- Radius Rescaling ----------------- //
   // ---------------------------------------------------- //
   
   if (rescale) {
  
     double density = tracers_catalogue->numdensity();
    
     cout << endl;
     coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;
  
     double rescaling_time = omp_get_wtime();
     
     vector<bool> remove(nObjects(), false), bad_rescale(nObjects(), false);
  
     double min_X = tracers_catalogue->Min(Var::_X_), max_X = tracers_catalogue->Max(Var::_X_);
     double min_Y = tracers_catalogue->Min(Var::_Y_), max_Y = tracers_catalogue->Max(Var::_Y_);
     double min_Z = tracers_catalogue->Min(Var::_Z_), max_Z = tracers_catalogue->Max(Var::_Z_);
  
 #pragma omp parallel num_threads(omp_get_max_threads())
     {     
 #pragma omp for ordered schedule(dynamic)
       for (size_t j=0; j<nObjects(); j++) {
         vector<double> temp_val = {3.*radius(j), xx(j)-min_X, max_X-xx(j), yy(j)-min_Y, max_Y-yy(j), zz(j)-min_Z, max_Z-zz(j), delta_r[1]};
         double value = *min_element(temp_val.begin(), temp_val.end());
 #pragma omp ordered
         ChM.get_searching_region(value);
         vector<long> close = ChM.close_objects(coordinate(j));
         vector<double> radii;
         for (auto&& k : close) {
           double distance = Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
           if (distance < value) radii.emplace_back(distance);
         }
         if (radii.size() < 3) {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
  
         std::sort(radii.begin(), radii.end());
         vector<vector<double>> data(2, vector<double>(radii.size())); //distances, density contrast 
         for (size_t i=0; i<radii.size()-1; i++) {
           data[0][i] = (radii[i]+radii[i+1])*0.5;
           data[1][i] = (i+1)/(volume_sphere(data[0][i])*density);
         }
         data[0][radii.size()-1] = 2*data[0][radii.size()-2] - data[0][radii.size()-3];
         data[1][radii.size()-1] = radii.size()/(volume_sphere(data[0][radii.size()-1])*density);
         
         unsigned int N = std::round(radii.size()*0.25);
         bool expand = false;
  
         if (N>1) {
           while (!(data[1][N-2] < threshold && data[1][N-1] > threshold)) 
             if (N>2) 
               N--;
             else {
               expand = true;
               break;
             }
           if (expand) {
             N = std::round(radii.size()*0.25);
             while (!(data[1][N-2] < threshold && data[1][N-1] > threshold)) 
               if (N<radii.size()-2) N++;
               else break;
           }
         }
       
         double new_radius;
         if (N<=2 or N>=radii.size()-2) {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
         else
           new_radius = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");
  
         if (new_radius > 0) set_var(j, Var::_Radius_, new_radius);
         else {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
  
     // -----------------------
  
         N = std::round(radii.size()*0.25);
         expand = false;
  
         if (N>1) {
           while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
         if (N>2)
           N--;
         else {
           expand = true;
           break;
         }
           if (expand) {
             N = std::round(radii.size()*0.25);
             while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
               if (N<radii.size()-2)
                 N++;
               else
                 break;
           }
         }
       
         double new_generic;
         if (N<=2 or N>=radii.size()-2) {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
         else
       new_generic = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");
  
         if (new_generic > 0) set_var(j, Var::_Generic_, new_generic);
         else {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
       }
     }
  
     remove_objects(remove);
     coutCBL << "Removed voids:" << endl;
     cout << "\t Bad rescaled: " << count(bad_rescale.begin(), bad_rescale.end(), true) << endl;
     
     vector<bool> remove_outofrange(nObjects(), false);
     for (size_t ii = 0; ii<nObjects(); ii++) 
       if (radius(ii) < delta_r[0] || delta_r[1] < radius(ii)) 
     remove_outofrange[ii] = true;
       
     remove_objects(remove_outofrange);
     cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"]: " << count(remove_outofrange.begin(), remove_outofrange.end(), true) << endl;
  
     compute_densityContrast(tracers_catalogue, ChM, ratio);
  
     coutCBL << "Time spent by the rescaling procedure: " << omp_get_wtime()-rescaling_time << " seconds" << endl;
     
   }
   
   coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
     
   // ---------------------------------------------------- //
   // ------------------ Overlap Check ------------------- //
   // ---------------------------------------------------- //
   
   if (checkoverlap) {
     cout << endl;
     coutCBL << par::col_green << " --- Checking for overlapping voids --- \n" << par::col_default << endl;
  
     double ol_time = omp_get_wtime();
     vector<bool> remove(nObjects(), false);
     vector<double> criteriumOrder(nObjects());
     string criterium;
  
     if (ol_criterion == Var::_CentralDensity_) {
       criteriumOrder = var(Var::_CentralDensity_);
       criterium = "central density";
     }
     else if (ol_criterion == Var::_DensityContrast_) {
       criteriumOrder = var(Var::_DensityContrast_);
       criterium = "density contrast";
     }
  
     else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");
  
     std::vector<int> indices(nObjects(), 0);
     std::iota(indices.begin(), indices.end(), 0);
     std::sort(indices.begin(), indices.end(), [&](int A, int B) 
           -> bool {return (ol_criterion == Var::_CentralDensity_) ? criteriumOrder[A] < criteriumOrder[B] : criteriumOrder[A] > criteriumOrder[B];});
     
     Order(indices);
     coutCBL << "Catalogue ordered according to the " << criterium << endl << endl;
     coutCBL << "* * * Generating ChainMesh for void centres * * *" << endl;
  
     chainmesh::ChainMesh3D ChM_voids(Min(Var::_Radius_), var(Var::_X_), var(Var::_Y_), var(Var::_Z_), 2*Max(Var::_Radius_));
  
     for (size_t i=0; i<nObjects(); i++) {
       if (!remove[i]) {
         vector<long> close = ChM_voids.close_objects(coordinate(i));
         std::sort(close.begin(), close.end());
         for (auto && j : close) {
           if (!remove[j]) {
             double distance = Catalogue::distance(i, catalogue_object(j));
             if (distance < radius(i)+radius(j) && (long)i!=j) {
               if (ol_criterion == Var::_CentralDensity_) {
                 if (centralDensity(i) < centralDensity(j)) {
                   remove[j] = true;
                 }
                 else if (centralDensity(i) > centralDensity(j)) {
                   remove[i] = true;
                   break;
                 }
                 else if (centralDensity(i) == centralDensity(j)) {
                   if (densityContrast(i) < densityContrast(j)) {
                     remove[i] = true;
                     break;
                   }
                   else 
                     remove[j] = true;
                 }
               }
               else if (ol_criterion == Var::_DensityContrast_) {
                 if (densityContrast(i) < densityContrast(j)) {
                   remove[i] = true;
                   break;
                 }
                 else if (densityContrast(i) > densityContrast(j)) 
                   remove[j] = true;
                 else if (densityContrast(i) == densityContrast(j)) {
                   if (centralDensity(i) < centralDensity(j)) 
                     remove[j] = true;
                   else {
                     remove[i] = true;
                     break;
                   }
                 }
               }
             }
           }
         }
       }
     }
  
     remove_objects(remove);
  
     cout << endl;
     coutCBL << "Voids removed to avoid overlap: " << count(remove.begin(), remove.end(), true) << endl;
     coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
   }
   
   cout << endl;
   coutCBL << "Voids in the final Catalogue: " << nObjects() << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
   
   cout << endl;
   coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;
  
 }
  
 //********************************** VERSION 2.0 **********************************//
  
 void cbl::catalogue::Catalogue::clean_void_catalogue (const std::vector<double> par_numdensity, const bool initial_radius, const std::vector<double> delta_r, const double threshold, bool rescale, const std::shared_ptr<Catalogue> tracers_catalogue, chainmesh::ChainMesh3D ChM, const double ratio, const bool checkoverlap, const Var ol_criterion)
 {
   const double start_time = omp_get_wtime();
     
   // ---------------------------------------------------- //
   // ---------------- Cleaning Procedure ---------------- //
   // ---------------------------------------------------- //
  
   cout << endl;
   coutCBL << par::col_blue << "> > > > > >>>>>>>>>> CLEANING PROCEDURE STARTED  <<<<<<<<< < < < < <" << par::col_default << endl << endl;
   coutCBL << "Voids in the initial Catalogue: " << nObjects() << endl << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
   coutCBL << par::col_green << " --- Removing spurious voids --- " << par::col_default << endl << endl;
   coutCBL << "Removed voids: " << endl;
   
   if (initial_radius) {
     vector<bool> remove(nObjects(), false);
     for (size_t i=0; i<nObjects(); i++) 
       if (radius(i) < delta_r[0] || delta_r[1] < radius(i)) remove[i] = true;
  
     remove_objects(remove);
     cout << "\t r_min-r_max criterion: " << count(remove.begin(), remove.end(), true) << endl;
   }
  
   compute_centralDensity(tracers_catalogue, ChM, par_numdensity, threshold);
   coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
  
   // ---------------------------------------------------- //
   // ----------------- Radius Rescaling ----------------- //
   // ---------------------------------------------------- //
   
   if (rescale) {
    
     cout << endl;
     coutCBL << par::col_green << " --- Rescaling radii --- \n" << par::col_default << endl;
  
     double rescaling_time = omp_get_wtime();
     
     vector<bool> remove(nObjects(), false), bad_rescale(nObjects(), false);
  
 #pragma omp parallel num_threads(omp_get_max_threads())
     {     
 #pragma omp for ordered schedule(dynamic)
       for (size_t j=0; j<nObjects(); j++) {
         double value = 2*radius(j);
 #pragma omp ordered
         ChM.get_searching_region(value);
         vector<long> close = ChM.close_objects(coordinate(j));
         vector<double> radii;
         for (auto&& k : close) {
           double distance = Catalogue::distance(j, tracers_catalogue->catalogue_object(k));
           if (distance < value) radii.emplace_back(distance);
         }
         if (radii.size() < 3) {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
  
         double zz = redshift(j);
         double density = 0.;
         for (size_t N=par_numdensity.size(); N-->0;) density += par_numdensity[par_numdensity.size()-1-N]*pow(zz,N);
  
         std::sort(radii.begin(), radii.end());
  
         vector<vector<double>> data(2, vector<double>(radii.size())); //distances, density contrast 
         for (size_t i=0; i<radii.size()-1; i++) {
           data[0][i] = (radii[i]+radii[i+1])*0.5;
           data[1][i] = (i+1)/(volume_sphere(data[0][i])*density);
         }
         data[0][radii.size()-1] = 2*data[0][radii.size()-2] - data[0][radii.size()-3];
         data[1][radii.size()-1] = radii.size()/(volume_sphere(data[0][radii.size()-1])*density);
         
         unsigned int N = std::round(radii.size()*0.25);
         bool expand = false;
  
         if (N>1) {
           while (!(data[1][N-2] < threshold && data[1][N-1] > threshold))
             if (N>2) N--;
             else {
               expand = true;
               break;
             }
           if (expand) {
             N = std::round(radii.size()*0.25);
             while (!(data[1][N-2] < threshold && data[1][N-1] > threshold)) 
               if (N<radii.size()-2) N++;
               else break;
           }
         }
       
         double new_radius;
         if (N<=2 or N>=radii.size()-2) {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
         else
           new_radius = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");
  
         if (new_radius > 0) set_var(j, Var::_Radius_, new_radius);
         else {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
  
     // -----------------------
  
         N = std::round(radii.size()*0.25);
         expand = false;
  
         if (N>1) {
           while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
         if (N>2)
           N--;
         else {
           expand = true;
           break;
         }
           if (expand) {
             N = std::round(radii.size()*0.25);
             while (!(data[1][N-2] < 0.5 && data[1][N-1] > 0.5))
               if (N<radii.size()-2)
                 N++;
               else
                 break;
           }
         }
       
         double new_generic;
         if (N<=2 or N>=radii.size()-2) {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
         else
       new_generic = interpolated(threshold, {data[1][N-1], data[1][N]}, {data[0][N-1], data[0][N]}, "Linear");
  
         if (new_generic > 0) set_var(j, Var::_Generic_, new_generic);
         else {
           remove[j] = true;
           bad_rescale[j] = true;
           continue;
         }
       }
     }
  
     remove_objects(remove);
     coutCBL << "Removed voids:" << endl;
     cout << "\t Bad rescaled: " << count(bad_rescale.begin(), bad_rescale.end(), true) << endl;
  
     vector<bool> remove_outofrange(nObjects(), false);
     for (size_t ii = 0; ii<nObjects(); ii++) 
       if (radius(ii) < delta_r[0] || delta_r[1] < radius(ii)) 
     remove_outofrange[ii] = true;
       
     remove_objects(remove_outofrange);
     cout << "\t Out of range ["+conv(delta_r[0],par::fDP2)+","+conv(delta_r[1],par::fDP2)+"]: " << count(remove_outofrange.begin(), remove_outofrange.end(), true) << endl;
     
     compute_densityContrast(tracers_catalogue, ChM, ratio);
  
     coutCBL << "Time spent by the rescaling procedure: " << omp_get_wtime()-rescaling_time << " seconds \n" << endl;
     
   }
   
   coutCBL << "Voids in the Catalogue: " << nObjects() << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp"); 
     
   // ---------------------------------------------------- //
   // ------------------ Overlap Check ------------------- //
   // ---------------------------------------------------- //
   
   if (checkoverlap) {
     cout << endl;
     coutCBL << par::col_green << " --- Checking for overlapping voids --- \n" << par::col_default << endl;
  
     double ol_time = omp_get_wtime();
     vector<bool> remove(nObjects(), false);
  
     vector<double> criteriumOrder(nObjects());
     string criterium;
  
     if (ol_criterion == Var::_CentralDensity_) {
       criteriumOrder = var(Var::_CentralDensity_);
       criterium = "central density";
     }
     else if (ol_criterion == Var::_DensityContrast_) {
       criteriumOrder = var(Var::_DensityContrast_);
       criterium = "density contrast";
     }
  
     else ErrorCBL("allowed overlap criteria are '_CentralDensity_' or '_DensityContrast_' .", "Catalogue", "VoidCatalogue.cpp");
  
     std::vector<int> indices(nObjects(), 0);
     std::iota(indices.begin(), indices.end(), 0);
     std::sort(indices.begin(), indices.end(), [&](int A, int B) 
           -> bool {return (ol_criterion == Var::_CentralDensity_) ? criteriumOrder[A] < criteriumOrder[B] : criteriumOrder[A] > criteriumOrder[B];});
     
     Order(indices);
     coutCBL << "Catalogue ordered according to the " << criterium << endl << endl;
     coutCBL << "* * * Generating ChainMesh for void centres * * *" << endl;
  
     chainmesh::ChainMesh3D ChM_voids(Min(Var::_Radius_), var(Var::_X_), var(Var::_Y_), var(Var::_Z_), 2*Max(Var::_Radius_));
     
     for (size_t i=0; i<nObjects(); i++) {
       if (!remove[i]) {
         vector<long> close = ChM_voids.close_objects(coordinate(i));
         std::sort(close.begin(), close.end());
         for (auto &&j : close) {
           if (!remove[j]) {
             double distance = Catalogue::distance(i, catalogue_object(j));
             if (distance<radius(i)+radius(j) && (long)i!=j) {
               if (ol_criterion == Var::_CentralDensity_) {
                 if (centralDensity(i) < centralDensity(j)) {
                   remove[j] = true;
                 }
                 else if (centralDensity(i) > centralDensity(j)) {
                   remove[i] = true;
                   break;
                 }
                 else if (centralDensity(i) == centralDensity(j)) {
                   if (densityContrast(i) < densityContrast(j)) {
                     remove[i] = true;
                     break;
                   }
                   else 
                     remove[j] = true;
                 }
               }
               else if (ol_criterion == Var::_DensityContrast_) {
                 if (densityContrast(i) < densityContrast(j)) {
                   remove[i] = true;
                   break;
                 }
                 else if (densityContrast(i) > densityContrast(j)) 
                   remove[j] = true;
                 else if (densityContrast(i) == densityContrast(j)) {
                   if (centralDensity(i) < centralDensity(j)) 
                     remove[j] = true;
                   else {
                     remove[i] = true;
                     break;
                   }
                 }
               }
             }
           }
         }
       }
     }
     
     remove_objects(remove);
  
     cout << endl;
     coutCBL << "Voids removed to avoid overlap: " << count(remove.begin(), remove.end(), true) << endl;
     coutCBL << "Time spent by the overlap-checking procedure: " << omp_get_wtime()-ol_time << " seconds" << endl;
   }
   
   cout << endl;
   coutCBL << "Voids in the final Catalogue: " << nObjects() << endl;
   if (nObjects()==0) ErrorCBL("Empty void catalogue!", "Catalogue", "VoidCatalogue.cpp");
   
   cout << endl;
   coutCBL << "Total time spent: " << omp_get_wtime()-start_time << " seconds \n" << endl;
  
 }