function [T_cond_ensemble, T_source_ensemble] = Tsite_Tsource_reconstruction_2020_ensemble(ice_core_dD, ice_core_d18O)
% cd '/Users/bradleymarkle/Documents/MATLAB/simple_water_isotope_model/'
cd '/Users/bradley/Documents/MATLAB/simple_water_isotope_model/'

method=1;
cd ./SWIM_results/2020/

dinfo = dir( fullfile('SWIM_results_2020*.mat') );
num_files = length(dinfo);
% filenames = fullfile( projectdir, {dinfo.name} );
filenames = fullfile({'./2020/'},{dinfo.name} );
% cd '/Users/bradleymarkle/Documents/MATLAB/simple_water_isotope_model/'
cd '/Users/bradley/Documents/MATLAB/simple_water_isotope_model'

for K = 1 : num_files
  this_file = filenames{K};
[T_cond_ensemble(:,K), T_source_ensemble(:,K)] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, this_file);

end

