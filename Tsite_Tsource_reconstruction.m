function [T_cond_ice_core, T_source_ice_core, T_surf_ice_core] = Tsite_Tsource_reconstruction(ice_core_dD, ice_core_d18O, method, file, run_model)
% [T_cond_ice_core, T_source_ice_core, T_surf_ice_core] = Tsite_Tsource_reconstruction(ice_core_dD, ice_core_d18O, method, file, run_model);
% This is a function that uses a simple water isotope model and a non-linear
% inversion to calculate T_source and T_site (both condensation temp and
% surface temp) for an ice core given the d18O and dD data.
% There are three (or more) methods that can be used. 
% Method = 1, it uses d18O_ln and dln for the inversion. (probably best)
% Method = 2, it uses d18O_ln and dD_ln.
% Method = 3, it uses d18O and dxs.
% If you want the program to run the model, set run_model = 1, if not let
% run_model = 0; The default is to not run the model.
% If you want to speficify a file of model data to use for the inversion,
% put it in input variable "file", e.g. file = 'SWIM_results/SWIM_results_test.mat';
% Otherwise the program will use a default file.
%% Housekeeping

if ~exist('method','var')
    method = 1;%default method is to use d18O and dxs
end

if ~exist('file','var')
    file = 'SWIM_results_local_ncep.mat';%default reconsturdtion file
%     file = 'SWIM_results_global_ncep.mat';%
%     file = 'SWIM_results_local_era.mat';%

end

if ~exist('run_model','var')
    run_model = 0;%default is to NOT run the model and use above file
end


%% set up ice core variables
ice_core_d18O_ln=(log(1 + ice_core_d18O./1000).*1000);
ice_core_dD_ln=(log(1 + ice_core_dD./1000).*1000);
ice_core_dln=ice_core_dD_ln-(-2.85*10^(-2) .* ice_core_d18O_ln.^2 + 8.47 .* ice_core_d18O_ln);
ice_core_dxs=ice_core_dD-8*ice_core_d18O;



%% either load data or run model
%% load existing run
if run_model ~= 1
    folder='./SWIM_results/';
    load(fullfile(folder,file))

% load('SWIM_results/SWIM_results_test2.mat')

%% run the model
elseif run_model == 1



Tsite_hi=-10; 
Tsite_lo=-51; 
dTsite=0.1; 
Tsite=[Tsite_lo Tsite_hi dTsite];

Tsource_hi=20; 
Tsource_lo=1; 
dTsource=0.1;
Tsource=[Tsource_lo Tsource_hi dTsource];

RHsource=[];
% RHsource=[0.1 0.9 0.1];

SH=1;
% closure = 'global';
closure = 'local';
reanalysis = 'ncep';
% reanalysis = 'era';
season = 'annual';
% season = 'DJF';
% season = 'JJA';


%tun the super saturation parameters
a=1;
b=0.00525;
% b=0.007;
% c=0.0000;
c=0.00001;


% [T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site] = simple_water_istope_model(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);
[T_site, T_source, RH_source, d18O_site, dD_site, d18Oln_site, dDln_site, dxs_site, d17O_xs_site, dlnU_site, r_s_site] = simple_water_istope_model_2018(Tsite, Tsource, RHsource, a, b, c, closure, reanalysis, SH, season);

% save('SWIM_results_test.mat','T_site','RH_source', 'T_source', 'd18O_site', 'dD_site', 'd18Oln_site', 'dDln_site', 'dxs_site', 'd17O_xs_site', 'dlnU_site');
end

%%

zz=repmat(T_source, length(T_site),1)';
% soursetest = griddata(d18Oln_site,dlnU_site,zz,-40,0)
ZZ=repmat(T_site, length(T_source),1);
% sitetest = griddata(d18Oln_site,dlnU_site,ZZ,-40,0)

if method == 1
T_cond_ice_core=griddata(d18Oln_site,dlnU_site,ZZ,ice_core_d18O_ln,ice_core_dln);
T_source_ice_core=griddata(d18Oln_site,dlnU_site,zz,ice_core_d18O_ln,ice_core_dln);
[T_surf_ice_core] = Ts_to_Tc(T_cond_ice_core,1);

elseif method == 2
T_cond_ice_core=griddata(d18Oln_site,dDln_site,ZZ,ice_core_d18O_ln,ice_core_dD_ln);
T_source_ice_core=griddata(d18Oln_site,dDln_site,zz,ice_core_d18O_ln,ice_core_dD_ln);
[T_surf_ice_core] = Ts_to_Tc(T_cond_ice_core,1);

elseif method == 3
T_cond_ice_core=griddata(d18Oln_site,dxs_site,ZZ,ice_core_d18O,ice_core_dxs);
T_source_ice_core=griddata(d18Oln_site,dxs_site,zz,ice_core_d18O,ice_core_dxs);
[T_surf_ice_core] = Ts_to_Tc(T_cond_ice_core,1);

end  
end

       