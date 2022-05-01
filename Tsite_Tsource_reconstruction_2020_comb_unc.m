function [T_cond_base, diff_total_cond, diff_total_cond_all, reldiff_total_cond, reldiff_total_cond_all, T_source_base, diff_total_source, diff_total_source_all, reldiff_total_source, reldiff_total_source_all, RH_offsets] = Tsite_Tsource_reconstruction_2020_comb_unc(ice_core_dD, ice_core_d18O)
method =1;


%% Set up files and reconstrucitons
% cd     '/Users/bradleymarkle/Documents/MATLAB/simple_water_isotope_model'
cd     '/Users/bradley/Documents/MATLAB/simple_water_isotope_model'

%set up base case
file_base= './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_base= './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_base= './2019/SWIM_results_20190828_ncep_annual_adiabat_local_initRH90_a1_b00525_c0_adiffevap_009.mat';
% file_base= './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initRH_a1_b00525_c0_adiffevap_009.mat';
% file_base= './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009.mat';
% file_base= './2019/SWIM_results_20190815_era_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009.mat';

[T_cond_base, T_source_base] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_base);

%%
%~~~~~~~~ a diff ~~~~~~~~
% file_adiff_007= './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_007.mat';
% file_adiff_008= './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_008.mat';
% file_adiff_009= './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009.mat';
% file_adiff_010= './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_010.mat';
file_adiff_007= './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_007_SOURCE_0_28_SITE_m70_0.mat';
file_adiff_008= './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_008_SOURCE_0_28_SITE_m70_0.mat';
file_adiff_009= './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_adiff_010= './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_010_SOURCE_0_28_SITE_m70_0.mat';

% [T_cond_adiff_007, T_source_adiff_007] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_adiff_007);
[T_cond_adiff_008, T_source_adiff_008] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_adiff_008);
[T_cond_adiff_009, T_source_adiff_009] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_adiff_009);
[T_cond_adiff_010, T_source_adiff_010] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_adiff_010);

T_cond_diff_adiff_1=T_cond_adiff_009-T_cond_adiff_008;
T_cond_diff_adiff_2=T_cond_adiff_009-T_cond_adiff_010;
T_cond_diff_adiff=nanmean([abs(T_cond_diff_adiff_1) abs(T_cond_diff_adiff_2)],2);

T_source_diff_adiff_1=T_source_adiff_009-T_source_adiff_008;
T_source_diff_adiff_2=T_source_adiff_009-T_source_adiff_010;
T_source_diff_adiff=nanmean([abs(T_source_diff_adiff_1) abs(T_source_diff_adiff_2)],2);


T_cond_reldiff_adiff_1=(T_cond_adiff_009-nanmean(T_cond_adiff_009))-(T_cond_adiff_008-nanmean(T_cond_adiff_008));
T_cond_reldiff_adiff_2=(T_cond_adiff_009-nanmean(T_cond_adiff_009))-(T_cond_adiff_010-nanmean(T_cond_adiff_010));
T_cond_reldiff_adiff=nanmean([abs(T_cond_reldiff_adiff_1) abs(T_cond_reldiff_adiff_2)],2);
T_source_reldiff_adiff_1=(T_source_adiff_009-nanmean(T_source_adiff_009))-(T_source_adiff_008-nanmean(T_source_adiff_008));
T_source_reldiff_adiff_2=(T_source_adiff_009-nanmean(T_source_adiff_009))-(T_source_adiff_010-nanmean(T_source_adiff_010));
T_source_reldiff_adiff=nanmean([abs(T_source_reldiff_adiff_1) abs(T_source_reldiff_adiff_2)],2);

%%
%~~~~~~~ tuning, b ~~~~~~~~~~~

% file_b00525 = './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009.mat';
% 
% % file_b005 = './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b005_c0_adiffevap_009.mat';
% file_b0051 = './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b0051_c0_adiffevap_009.mat';
% file_b0052 = './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b0052_c0_adiffevap_009.mat';
% file_b0053 = './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b0053_c0_adiffevap_009.mat';
% file_b0054 = './2019/SWIM_results_20190514_ncep_annual_adiabat_local_initsat_a1_b0054_c0_adiffevap_009.mat';
% 
% file_b00525_c = './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initsat_a1_b00525_c00001_adiffevap_009.mat';

file_b00525 = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_b0051 = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b0051_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_b0052 = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b0052_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_b0053 = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b0053_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_b0054 = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b0054_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_b00525_c = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c00001_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';


% file_b00525 = './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_b0051 = './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00510_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_b0052 = './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00520_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_b0053 = './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00530_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_b0054 = './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00540_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';

[T_cond_b00525, T_source_b00525] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_b00525);
% [T_cond_b005, T_source_b005] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_b005);
[T_cond_b0051, T_source_b0051] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_b0051);
[T_cond_b0052, T_source_b0052] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_b0052);
[T_cond_b0053, T_source_b0053] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_b0053);
[T_cond_b0054, T_source_b0054] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_b0054);
[T_cond_b00525_c, T_source_b00525_c] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_b00525_c);

T_cond_diff_b_1=T_cond_b00525-T_cond_b0051;
T_cond_diff_b_2=T_cond_b00525-T_cond_b0052;
T_cond_diff_b_3=T_cond_b00525-T_cond_b0053;
T_cond_diff_b_4=T_cond_b00525-T_cond_b0054;
% T_cond_diff_b_c=T_cond_b00525-T_cond_b00525_c;
% T_cond_diff_b=nanmean([abs(T_cond_diff_b_1) abs(T_cond_diff_b_2) abs(T_cond_diff_b_3) abs(T_cond_diff_b_4) abs(T_cond_diff_b_c)]);
T_cond_diff_b=nanmean([abs(T_cond_diff_b_1) abs(T_cond_diff_b_2) abs(T_cond_diff_b_3) abs(T_cond_diff_b_4)],2);

T_source_diff_b_1=T_source_b00525-T_source_b0051;
T_source_diff_b_2=T_source_b00525-T_source_b0052;
T_source_diff_b_3=T_source_b00525-T_source_b0053;
T_source_diff_b_4=T_source_b00525-T_source_b0054;
% T_source_diff_b_c=T_source_b00525-T_source_b00525_c;
% T_source_diff_b=nanmean([abs(T_source_diff_b_1) abs(T_source_diff_b_2) abs(T_source_diff_b_3) abs(T_source_diff_b_4) abs(T_source_diff_b_c)]);
T_source_diff_b=nanmean([abs(T_source_diff_b_1) abs(T_source_diff_b_2) abs(T_source_diff_b_3) abs(T_source_diff_b_4)],2);


T_cond_reldiff_b_1=(T_cond_b00525-nanmean(T_cond_b00525))-(T_cond_b0051-nanmean(T_cond_b0051));
T_cond_reldiff_b_2=(T_cond_b00525-nanmean(T_cond_b00525))-(T_cond_b0052-nanmean(T_cond_b0052));
T_cond_reldiff_b_3=(T_cond_b00525-nanmean(T_cond_b00525))-(T_cond_b0053-nanmean(T_cond_b0053));
T_cond_reldiff_b_4=(T_cond_b00525-nanmean(T_cond_b00525))-(T_cond_b0054-nanmean(T_cond_b0054));
% T_cond_reldiff_b_c=(T_cond_b00525-nanmean(T_cond_b00525))-(T_cond_b0525_c-nanmean(T_cond_b0525_c));
% T_cond_reldiff_b=nanmean([abs(T_cond_reldiff_b_1) abs(T_cond_reldiff_b_2) abs(T_cond_reldiff_b_3) abs(T_cond_reldiff_b_4) abs(T_cond_reldiff_b_c)]);
T_cond_reldiff_b=nanmean([abs(T_cond_reldiff_b_1) abs(T_cond_reldiff_b_2) abs(T_cond_reldiff_b_3) abs(T_cond_reldiff_b_4)],2);

T_source_reldiff_b_1=(T_source_b00525-nanmean(T_source_b00525))-(T_source_b0051-nanmean(T_source_b0051));
T_source_reldiff_b_2=(T_source_b00525-nanmean(T_source_b00525))-(T_source_b0052-nanmean(T_source_b0052));
T_source_reldiff_b_3=(T_source_b00525-nanmean(T_source_b00525))-(T_source_b0053-nanmean(T_source_b0053));
T_source_reldiff_b_4=(T_source_b00525-nanmean(T_source_b00525))-(T_source_b0054-nanmean(T_source_b0054));
% T_source_reldiff_b_c=(T_source_b00525-nanmean(T_source_b00525))-(T_source_b0525_c-nanmean(T_source_b0525_c));
% T_source_reldiff_b=nanmean([abs(T_source_reldiff_b_1) abs(T_source_reldiff_b_2) abs(T_source_reldiff_b_3) abs(T_source_reldiff_b_4) abs(T_source_reldiff_b_c)]);
T_source_reldiff_b=nanmean([abs(T_source_reldiff_b_1) abs(T_source_reldiff_b_2) abs(T_source_reldiff_b_3) abs(T_source_reldiff_b_4)],2);

% figure;
% hold on
% plot(T_source_diff_b)
% plot(T_source_reldiff_b)

%%
%~~~~~~~ closure assumption and reanalysis

% file_ncep_global = './2019/SWIM_results_20190815_ncep_annual_adiabat_global_initsat_a1_b00525_c0_adiffevap_009.mat';
% file_ncep_local = './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009.mat';
% % file_era_local = './2019/SWIM_results_20190514_era_annual_adiabat_local_initRH_a1_b00525_c0_adiffevap_009.mat';
% % file_era_global = './2019/SWIM_results_20190514_era_annual_adiabat_global_initRH_a1_b00525_c0_adiffevap_009.mat';
% file_era_local2 = './2019/SWIM_results_20190815_era_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009.mat';
% file_era_global2= './2019/SWIM_results_20190514_era_annual_adiabat_global_initsat_a1_b00525_c0_adiffevap_009.mat';

file_ncep_global = './2020/SWIM_results_202006_ncep_annual_adiabat_global_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_ncep_local = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_era_local = './2019/SWIM_results_20190514_era_annual_adiabat_local_initRH_a1_b00525_c0_adiffevap_009.mat';
% file_era_global = './2019/SWIM_results_20190514_era_annual_adiabat_global_initRH_a1_b00525_c0_adiffevap_009.mat';
file_era_local2 = './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_era_global2= './2020/SWIM_results_202006_era_annual_adiabat_global_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';


[T_cond_ncep_global, T_source_ncep_global] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_ncep_global);
[T_cond_ncep_local, T_source_ncep_local] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_ncep_local);
% [T_cond_era_global, T_source_era_global] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_era_global);
% [T_cond_era_local, T_source_era_local] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_era_local);
[T_cond_era_global2, T_source_era_global2] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_era_global2);
[T_cond_era_local2, T_source_era_local2] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_era_local2);

% figure;hold on
% plot(T_source_ncep_global)
% plot(T_source_ncep_local)
% plot(T_source_era_global)
% plot(T_source_era_local)
% plot(T_source_era_global2)
% plot(T_source_era_local2)

%closure
T_cond_diff_closure_1=T_cond_ncep_local-T_cond_ncep_global;
% T_cond_diff_closure_2=T_cond_era_local-T_cond_era_global;
T_cond_diff_closure_3=T_cond_era_local2-T_cond_era_global2;
T_cond_diff_closure=nanmean([abs(T_cond_diff_closure_1) abs(T_cond_diff_closure_3)],2);

T_source_diff_closure_1=T_source_ncep_local-T_source_ncep_global;
% T_source_diff_closure_2=T_source_era_local-T_source_era_global;
T_source_diff_closure_3=T_source_era_local2-T_source_era_global2;
T_source_diff_closure=nanmean([abs(T_source_diff_closure_1) abs(T_source_diff_closure_3)],2);

T_cond_reldiff_closure_1=(T_cond_ncep_local-nanmean(T_cond_ncep_local))-(T_cond_ncep_global-nanmean(T_cond_ncep_global));
% T_cond_reldiff_closure_2=(T_cond_era_local-nanmean(T_cond_era_local))-(T_cond_era_global-nanmean(T_cond_era_global));
T_cond_reldiff_closure_3=(T_cond_era_local2-nanmean(T_cond_era_local2))-(T_cond_era_global2-nanmean(T_cond_era_global2));
T_cond_reldiff_closure=nanmean([abs(T_cond_reldiff_closure_1) abs(T_cond_reldiff_closure_3)],2);

T_source_reldiff_closure_1=(T_source_ncep_local-nanmean(T_source_ncep_local))-(T_source_ncep_global-nanmean(T_source_ncep_global));
% T_source_reldiff_closure_2=(T_source_era_local-nanmean(T_source_era_local))-(T_source_era_global-nanmean(T_source_era_global));
T_source_reldiff_closure_3=(T_source_era_local2-nanmean(T_source_era_local2))-(T_source_era_global2-nanmean(T_source_era_global2));
T_source_reldiff_closure=nanmean([abs(T_source_reldiff_closure_1) abs(T_source_reldiff_closure_3)],2);


%reanalysis
T_cond_diff_reanalysis_1=T_cond_ncep_local-T_cond_era_local2;
T_cond_diff_reanalysis_2=T_cond_ncep_global-T_cond_era_global2;
T_cond_diff_reanalysis=nanmean([abs(T_cond_diff_reanalysis_1) abs(T_cond_diff_reanalysis_2)],2);

T_source_diff_reanalysis_1=T_source_ncep_local-T_source_era_local2;
T_source_diff_reanalysis_2=T_source_ncep_global-T_source_era_global2;
T_source_diff_reanalysis=nanmean([abs(T_source_diff_reanalysis_1) abs(T_source_diff_reanalysis_2)],2);

T_cond_reldiff_reanalysis_1=(T_cond_ncep_local-nanmean(T_cond_ncep_local))-(T_cond_era_local2-nanmean(T_cond_era_local2));
T_cond_reldiff_reanalysis_2=(T_cond_ncep_global-nanmean(T_cond_ncep_global))-(T_cond_era_global2-nanmean(T_cond_era_global2));
T_cond_reldiff_reanalysis=nanmean([abs(T_cond_reldiff_reanalysis_1) abs(T_cond_reldiff_reanalysis_2)],2);

T_source_reldiff_reanalysis_1=(T_source_ncep_local-nanmean(T_source_ncep_local))-(T_source_era_local2-nanmean(T_source_era_local2));
T_source_reldiff_reanalysis_2=(T_source_ncep_global-nanmean(T_source_ncep_global))-(T_source_era_global2-nanmean(T_source_era_global2));
T_source_reldiff_reanalysis=nanmean([abs(T_source_reldiff_reanalysis_1) abs(T_source_reldiff_reanalysis_2)],2);

%%
%~~~~~~~
%~~~~~~~ precipitaion parameterization
% file_initRH80 = './2019/SWIM_results_20190828_ncep_annual_adiabat_local_initRH80_a1_b00525_c0_adiffevap_009.mat';
% file_initRH90 = './2019/SWIM_results_20190828_ncep_annual_adiabat_local_initRH90_a1_b00525_c0_adiffevap_009.mat';
% file_initRH = './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initRH_a1_b00525_c0_adiffevap_009.mat';
% file_initsat = './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009.mat';
file_initRH80 = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initRH80_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_initRH90 = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initRH90_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_initRH = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initRH_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
file_initsat = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';




[T_cond_initRH80, T_source_initRH80] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_initRH80);
[T_cond_initRH90, T_source_initRH90] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_initRH90);
[T_cond_initRH, T_source_initRH] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_initRH);
[T_cond_initsat, T_source_initsat] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_initsat);

T_cond_diff_precip1=T_cond_initRH90-T_cond_initRH80;
T_cond_diff_precip2=T_cond_initRH90-T_cond_initRH;
T_cond_diff_precip3=T_cond_initRH90-T_cond_initsat;
T_cond_diff_precip=nanmean([abs(T_cond_diff_precip1) abs(T_cond_diff_precip2) abs(T_cond_diff_precip3)],2);

T_source_diff_precip1=T_source_initRH90-T_source_initRH80;
T_source_diff_precip2=T_source_initRH90-T_source_initRH;
T_source_diff_precip3=T_source_initRH90-T_source_initsat;
T_source_diff_precip=nanmean([abs(T_source_diff_precip1) abs(T_source_diff_precip2) abs(T_source_diff_precip3)],2);

T_cond_reldiff_precip1=(T_cond_initRH90-nanmean(T_cond_initRH90))-(T_cond_initRH80-nanmean(T_cond_initRH80));
T_cond_reldiff_precip2=(T_cond_initRH90-nanmean(T_cond_initRH90))-(T_cond_initRH-nanmean(T_cond_initRH));
T_cond_reldiff_precip3=(T_cond_initRH90-nanmean(T_cond_initRH90))-(T_cond_initsat-nanmean(T_cond_initsat));
T_cond_reldiff_precip=nanmean([abs(T_cond_reldiff_precip1) abs(T_cond_reldiff_precip2) abs(T_cond_reldiff_precip3)],2);

T_source_reldiff_precip1=(T_source_initRH90-nanmean(T_source_initRH90))-(T_source_initRH80-nanmean(T_source_initRH80));
T_source_reldiff_precip2=(T_source_initRH90-nanmean(T_source_initRH90))-(T_source_initRH-nanmean(T_source_initRH));
T_source_reldiff_precip3=(T_source_initRH90-nanmean(T_source_initRH90))-(T_source_initsat-nanmean(T_source_initsat));
T_source_reldiff_precip=nanmean([abs(T_source_reldiff_precip1) abs(T_source_reldiff_precip2) abs(T_source_reldiff_precip3)],2);

%%
%~~~~~~~
%~~~~~~~ initial relative humidity
RH_offsets=0;
%off for now but look at later.
% 
% file_RH_80 = './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_RH80.mat';
% file_RH_85 = './2019/SWIM_results_20190815_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_RH85.mat';
% file_RH_75 = './2019/SWIM_results_20190904_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_RH75.mat';
% 
% [T_cond_RH_75, T_source_RH_75] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_RH_75);
% [T_cond_RH_80, T_source_RH_80] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_RH_80);
% [T_cond_RH_85, T_source_RH_85] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_RH_85);
% 
% T_cond_diff_RH_1=T_cond_initsat-T_cond_RH_75;
% T_cond_diff_RH_2=T_cond_initsat-T_cond_RH_80;
% T_cond_diff_RH_3=T_cond_initsat-T_cond_RH_85;
% T_source_diff_RH_1=T_source_initsat-T_source_RH_75;
% T_source_diff_RH_2=T_source_initsat-T_source_RH_80;
% T_source_diff_RH_3=T_source_initsat-T_source_RH_85;
% 
% % RH_offsets = [T_cond_diff_RH_1 T_cond_diff_RH_2 T_cond_diff_RH_3 T_source_diff_RH_1 T_source_diff_RH_2 T_source_diff_RH_3];
% 
% file_initRH90_minus5 = './2019/SWIM_results_20190904_ncep_annual_adiabat_local_initRH90_a1_b00525_c0_adiffevap_009_RHoffsetminus5.mat';
% file_initRH90_minus10 = './2019/SWIM_results_20190904_ncep_annual_adiabat_local_initRH90_a1_b00525_c0_adiffevap_009_RHoffsetminus10.mat';
% [T_cond_initRH90_minus5, T_source_initRH90_minus5] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_initRH90_minus5);
% [T_cond_initRH90_minus10, T_source_initRH90_minus10] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_initRH90_minus10);
% 
% T_cond_diff_RH_alt=T_cond_initRH90-T_cond_initRH90_minus5;
% T_cond_diff_RH_alt2=T_cond_initRH90-T_cond_initRH90_minus10;
% 
% T_source_diff_RH_alt=T_source_initRH90-T_source_initRH90_minus5;
% T_source_diff_RH_alt2=T_source_initRH90-T_source_initRH90_minus10;
% 
% RH_offsets = [T_cond_diff_RH_1 T_cond_diff_RH_2 T_cond_diff_RH_3 T_cond_diff_RH_alt T_cond_diff_RH_alt2...
%     T_source_diff_RH_1 T_source_diff_RH_2 T_source_diff_RH_3 T_source_diff_RH_alt T_source_diff_RH_alt2];

%~~~~~~~
%~~~~~~~


% 
% figure;hold on
% plot(T_cond_initRH90)
% plot(T_cond_initRH90_minus5)
% plot(T_cond_initRH90_minus10)
% 
% figure;hold on
% plot(T_source_initRH90)
% plot(T_source_initRH90_minus5)
% plot(T_source_initRH90_minus10)
% 
% figure;hold on
% plot(T_source_diff_RH_1)
% plot(T_source_diff_RH_2)
% plot(T_source_diff_RH_3)
% plot(T_source_diff_RH_alt)
% plot(T_source_diff_RH_alt2)
% 
% 
% figure;hold on
% plot(T_cond_diff_RH_1)
% plot(T_cond_diff_RH_2)
% plot(T_cond_diff_RH_3)
% plot(T_cond_diff_RH_alt)
% plot(T_cond_diff_RH_alt2)
% 
% 
% % 
% figure;hold on;
% % plot(T_source_RH_75);
% % plot(T_source_RH_80);
% % plot(T_source_RH_85);
% plot(T_source_initRH)
% plot(T_source_initRH80,'c')
% plot(T_source_initRH90,'b')
% plot(T_source_initsat,'k')
% plot(T_source_initRH90)
% plot(T_source_initRH90_minus5)
% plot(T_source_initRH90_minus10)
% 
% 
% figure;hold on;
% % xxxx=[T_source_initRH90 T_source_initRH90_minus5 T_source_initsat T_source_initRH90 T_source_initRH80 T_source_initRH];
% xxxx=[T_source_initRH90 T_source_initRH90_minus5 T_source_initRH90_minus10 T_source_initsat T_source_initRH90 T_source_initRH80 T_source_initRH T_source_RH_75 T_source_RH_80 T_source_RH_85];
% plot(mean(xxxx,2)-std(xxxx,[],2),'c')
% plot(mean(xxxx,2)+std(xxxx,[],2),'c')
% plot(mean(xxxx,2))


%~~~~~~~
%~~~~~~~ spatial variability in seawater isotope values
% !!!!! This one is tricky due to temporal variation in d18O_SW!!! leave
% for later
% file_fixed_d18Osw = './2019/SWIM_results_20190828_ncep_annual_adiabat_local_initRH90_a1_b00525_c0_adiffevap_009.mat';
% file_var_d18Osw = './2019/SWIM_results_20190903_ncep_annual_adiabat_local_initRH90_a1_b00525_c0_adiffevap_009_d18Osw_minus03_03.mat';
% 
% 
% [T_cond_fixedSW, T_source_fixedSW] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_fixed_d18Osw);
% [T_cond_varSW, T_source_varSW] = Tsite_Tsource_reconstruction_quick(ice_core_dD, ice_core_d18O, method, file_var_d18Osw );
% 
% T_cond_diff_SW=T_cond_fixedSW-T_cond_varSW;
% T_source_diff_SW=T_source_fixedSW-T_source_varSW;


%% Mixing uncertainty

%cond diff 0.055
%souorce diff 0.0197
%cond reldiff 0.159
%sourcd reldfif 0.32
%% nonuniqueness

% T_source_diff_nonunique

%% combine all uncertainties in quadrature


% diff_total_source=((T_source_diff_closure./2).^2 + (T_source_diff_adiff).^2 +(T_source_diff_b).^2 +(T_source_diff_precip).^2 +(T_source_diff_reanalysis).^2  ).^(1/2);
% diff_total_cond=((T_cond_diff_closure./2).^2 + (T_cond_diff_adiff).^2 +(T_cond_diff_b).^2 +(T_cond_diff_precip).^2 +(T_cond_diff_reanalysis).^2  ).^(1/2);
% 
% reldiff_total_source=((T_source_reldiff_closure./2).^2 + (T_source_reldiff_adiff).^2 +(T_source_reldiff_b).^2 +(T_source_reldiff_precip).^2 +(T_source_reldiff_reanalysis).^2  ).^(1/2);
% reldiff_total_cond=((T_cond_reldiff_closure./2).^2 + (T_cond_reldiff_adiff).^2 +(T_cond_reldiff_b).^2 +(T_cond_reldiff_precip).^2 +(T_cond_reldiff_reanalysis).^2  ).^(1/2);
% 

% %~~~~~~ put all uncertainties together
% diff_total_cond_all=[(T_cond_diff_closure./2)  (T_cond_diff_adiff) (T_cond_diff_b) (T_cond_diff_precip) (T_cond_diff_reanalysis)];
% diff_total_source_all=[(T_source_diff_closure./2)  (T_source_diff_adiff) (T_source_diff_b) (T_source_diff_precip) (T_source_diff_reanalysis)];
% 
% reldiff_total_cond_all=[(T_cond_reldiff_closure./2)  (T_cond_reldiff_adiff) (T_cond_reldiff_b) (T_cond_reldiff_precip) (T_cond_reldiff_reanalysis)];
% reldiff_total_source_all=[(T_source_reldiff_closure./2)  (T_source_reldiff_adiff) (T_source_reldiff_b) (T_source_reldiff_precip) (T_source_reldiff_reanalysis)];

% %~~~~~~try only closure adiff and b
% diff_total_cond_all=[(T_cond_diff_closure./2)  (T_cond_diff_adiff) (T_cond_diff_b)];
% diff_total_source_all=[(T_source_diff_closure./2)  (T_source_diff_adiff) (T_source_diff_b) ];
% 
% reldiff_total_cond_all=[(T_cond_reldiff_closure./2)  (T_cond_reldiff_adiff) (T_cond_reldiff_b)];
% reldiff_total_source_all=[(T_source_reldiff_closure./2)  (T_source_reldiff_adiff) (T_source_reldiff_b)];

% %~~~~~~ put all uncertainties together, plus add relative uncertainty due to
% % mixing!! 
diff_total_cond_all=[(T_cond_diff_closure./2)  (T_cond_diff_adiff) (T_cond_diff_b) (T_cond_diff_precip) (T_cond_diff_reanalysis) 0.055.*ones(size(T_cond_diff_reanalysis))];
diff_total_source_all=[(T_source_diff_closure./2)  (T_source_diff_adiff) (T_source_diff_b) (T_source_diff_precip) (T_source_diff_reanalysis) 0.0197.*ones(size(T_source_diff_reanalysis))];

reldiff_total_cond_all=[(T_cond_reldiff_closure./2)  (T_cond_reldiff_adiff) (T_cond_reldiff_b) (T_cond_reldiff_precip) (T_cond_reldiff_reanalysis) 0.159.*ones(size(T_cond_reldiff_reanalysis))];
reldiff_total_source_all=[(T_source_reldiff_closure./2)  (T_source_reldiff_adiff) (T_source_reldiff_b) (T_source_reldiff_precip) (T_source_reldiff_reanalysis) 0.32.*ones(size(T_source_reldiff_reanalysis))];


diff_total_cond=sum(diff_total_cond_all.^2,2).^(1/2);
diff_total_source=sum(diff_total_source_all.^2,2).^(1/2);

reldiff_total_cond=sum(reldiff_total_cond_all.^2,2).^(1/2);
reldiff_total_source=sum(reldiff_total_source_all.^2,2).^(1/2);

% 
% figure;
% hold on
% plot(reldiff_total_source)
% plot(diff_total_source)

% figure;
% hold on
% plot(reldiff_total_cond)
% plot(diff_total_cond)
% 
% 
% figure;
% subplot(211)
% hold on
% plot(T_cond_base+diff_total_cond,'c')
% plot(T_cond_base-diff_total_cond,'c')
% plot(T_cond_base+reldiff_total_cond,'b')
% plot(T_cond_base-reldiff_total_cond,'b')
% plot(T_cond_base,'r')
% subplot(212)
% hold on
% plot(T_source_base+diff_total_source,'c')
% plot(T_source_base-diff_total_source,'c')
% plot(T_source_base+reldiff_total_source,'b')
% plot(T_source_base-reldiff_total_source,'b')
% plot(T_source_base,'r')


end
