%% Load Data
% cd('/Users/Bradley/Documents/work/Antarctic_temperature/ice_core_data')
% cd('/Volumes/Macintosh HD/Users/Bradley/Documents/work/Antarctic_temperature/ice_core_data')
cd '/Users/bradley/work/projects/Antarctic_temperature/ice_core_data'
% compile_deep_data_2018
compile_deep_data_2020_synched
% compile_deep_data_2020_unsynched

%% correct for seawater
cd('/Users/bradley/Documents/MATLAB/simple_water_isotope_model/')

d18Osw_i=-0.3;
[wais_d18O_ln_cor,wais_dD_ln_cor,wais_dln_cor, wais_d18O_cor, wais_dD_cor]= seawater_cor_ln(wais_d18O,wais_dD,wais_age',d18Osw_i);
[edc_d18O_ln_cor,edc_dD_ln_cor,edc_dln_cor, edc_d18O_cor, edc_dD_cor]= seawater_cor_ln(edc_d18O,edc_dD,edc_age',d18Osw_i);
[edml_d18O_ln_cor,edml_dD_ln_cor,edml_dln_cor, edml_d18O_cor, edml_dD_cor]= seawater_cor_ln(edml_d18O,edml_dD,edml_age',d18Osw_i);
[vostok_d18O_ln_cor,vostok_dD_ln_cor,vostok_dln_cor, vostok_d18O_cor, vostok_dD_cor]= seawater_cor_ln(vostok_d18O,vostok_dD,vostok_age,d18Osw_i);
[siple_d18O_ln_cor,siple_dD_ln_cor,siple_dln_cor, siple_d18O_cor, siple_dD_cor]= seawater_cor_ln(siple_d18O,siple_dD,siple_age,d18Osw_i);
[talos_d18O_ln_cor,talos_dD_ln_cor,talos_dln_cor, talos_d18O_cor, talos_dD_cor]= seawater_cor_ln(talos_d18O,talos_dD,talos_age',d18Osw_i);
[taylor_d18O_ln_cor,taylor_dD_ln_cor,taylor_dln_cor, taylor_d18O_cor, taylor_dD_cor]= seawater_cor_ln(taylor_d18O,taylor_dD,taylor_age,d18Osw_i);
[fuji_d18O_ln_cor,fuji_dD_ln_cor,fuji_dln_cor, fuji_d18O_cor, fuji_dD_cor]= seawater_cor_ln(fuji_d18O,fuji_dD,fuji_age',d18Osw_i);
% [fuji2_d18O_ln_cor,fuji2_dD_ln_cor,fuji2_dln_cor, fuji2_d18O_cor, fuji2_dD_cor]= seawater_cor_ln(fuji2_d18O,fuji2_dD,fuji2_age',d18Osw_i);

[SP_d18O_ln_cor,SP_dD_ln_cor,SP_dln_cor, SP_d18O_cor, SP_dD_cor]= seawater_cor_ln(SP_d18O,SP_dD,SP_age,d18Osw_i);



%%
cd('/Users/bradley/Documents/MATLAB/simple_water_isotope_model/SWIM_results')


method = 1;


file_base = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';

% file_base = './2020/SWIM_results_202006_era_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';

cd('/Users/bradley/Documents/MATLAB/simple_water_isotope_model/')


%%

%% Base Case
% [wais_T_site, wais_T_source, wais_T_surf] = Tsite_Tsource_reconstruction_2019(wais_dD, wais_d18O,method,file_base);
% [wais_T_cond_cor, wais_T_source_cor, wais_T_surf_cor] = Tsite_Tsource_reconstruction_2019(wais_dD_cor, wais_d18O_cor,method,file_base);
% 
% [edc_T_cond, edc_T_source, edc_T_surf] = Tsite_Tsource_reconstruction_2019(edc_dD, edc_d18O,method,file_base);
% [edc_T_cond_cor, edc_T_source_cor, edc_T_surf_cor] = Tsite_Tsource_reconstruction_2019(edc_dD_cor, edc_d18O_cor,method,file_base);
% 
% 
% [edml_T_cond, edml_T_source, edml_T_surf] = Tsite_Tsource_reconstruction_2019(edml_dD, edml_d18O,method,file_base);
% [edml_T_cond_cor, edml_T_source_cor, edml_T_surf_cor] = Tsite_Tsource_reconstruction_2019(edml_dD_cor, edml_d18O_cor,method,file_base);
% 
% [vostok_T_cond, vostok_T_source, vostok_T_surf] = Tsite_Tsource_reconstruction_2019(vostok_dD, vostok_d18O,method,file_base);
% [vostok_T_cond_cor, vostok_T_source_cor, vostok_T_surf_cor] = Tsite_Tsource_reconstruction_2019(vostok_dD_cor, vostok_d18O_cor,method,file_base);
% 
% [siple_T_cond, siple_T_source, siple_T_surf] = Tsite_Tsource_reconstruction_2019(siple_dD, siple_d18O,method,file_base);
% [siple_T_cond_cor, siple_T_source_cor, siple_T_surf_cor] = Tsite_Tsource_reconstruction_2019(siple_dD_cor, siple_d18O_cor,method,file_base);
% 
% [talos_T_cond, talos_T_source, talos_T_surf] = Tsite_Tsource_reconstruction_2019(talos_dD, talos_d18O,method,file_base);
% [talos_T_cond_cor, talos_T_source_cor, talos_T_surf_cor] = Tsite_Tsource_reconstruction_2019(talos_dD_cor, talos_d18O_cor,method,file_base);
% 
% [fuji_T_cond, fuji_T_source, fuji_T_surf] = Tsite_Tsource_reconstruction_2019(fuji_dD, fuji_d18O,method,file_base);
% [fuji_T_cond_cor, fuji_T_source_cor, fuji_T_surf_cor] = Tsite_Tsource_reconstruction_2019(fuji_dD_cor, fuji_d18O_cor,method,file_base);
% 
% % [fuji2_T_cond, fuji2_T_source, fuji2_T_surf] = Tsite_Tsource_reconstruction_2019(fuji2_dD, fuji2_d18O,method,file_base);
% % [fuji2_T_cond_cor, fuji2_T_source_cor, fuji2_T_surf_cor] = Tsite_Tsource_reconstruction_2019(fuji2_dD_cor, fuji2_d18O_cor,method,file_base);
% 
% [SP_T_cond, SP_T_source, SP_T_surf] = Tsite_Tsource_reconstruction_2019(SP_dD, SP_d18O,method,file_base);
% [SP_T_cond_cor, SP_T_source_cor, SP_T_surf_cor] = Tsite_Tsource_reconstruction_2019(SP_dD_cor, SP_d18O_cor,method,file_base);
% 
% [taylor_T_cond, taylor_T_source, taylor_T_surf] = Tsite_Tsource_reconstruction_2019(taylor_dD, taylor_d18O,method,file_base);
% [taylor_T_cond_cor, taylor_T_source_cor, taylor_T_surf_cor] = Tsite_Tsource_reconstruction_2019(taylor_dD_cor, taylor_d18O_cor,method,file_base);
% 
% [talos_T_cond, talos_T_source, talos_T_surf] = Tsite_Tsource_reconstruction_2019(talos_dD, talos_d18O,method,file_base);
% [talos_T_cond_cor, talos_T_source_cor, talos_T_surf_cor] = Tsite_Tsource_reconstruction_2019(talos_dD_cor, talos_d18O_cor,method,file_base);

[wais_T_cond, wais_T_source, wais_rs] = Tsite_Tsource_reconstruction_quick(wais_dD, wais_d18O,method,file_base);
[wais_T_cond_cor, wais_T_source_cor, wais_rs_cor] = Tsite_Tsource_reconstruction_quick(wais_dD_cor, wais_d18O_cor,method,file_base);

[edc_T_cond, edc_T_source, edc_rs] = Tsite_Tsource_reconstruction_quick(edc_dD, edc_d18O,method,file_base);
[edc_T_cond_cor, edc_T_source_cor, edc_rs_cor] = Tsite_Tsource_reconstruction_quick(edc_dD_cor, edc_d18O_cor,method,file_base);


[edml_T_cond, edml_T_source, edml_rs] = Tsite_Tsource_reconstruction_quick(edml_dD, edml_d18O,method,file_base);
[edml_T_cond_cor, edml_T_source_cor, edml_rs_cor] = Tsite_Tsource_reconstruction_quick(edml_dD_cor, edml_d18O_cor,method,file_base);

[vostok_T_cond, vostok_T_source, vostok_rs] = Tsite_Tsource_reconstruction_quick(vostok_dD, vostok_d18O,method,file_base);
[vostok_T_cond_cor, vostok_T_source_cor, vostok_rs_cor] = Tsite_Tsource_reconstruction_quick(vostok_dD_cor, vostok_d18O_cor,method,file_base);

[siple_T_cond, siple_T_source, siple_rs] = Tsite_Tsource_reconstruction_quick(siple_dD, siple_d18O,method,file_base);
[siple_T_cond_cor, siple_T_source_cor, siple_rs_cor] = Tsite_Tsource_reconstruction_quick(siple_dD_cor, siple_d18O_cor,method,file_base);

[fuji_T_cond, fuji_T_source, fuji_rs] = Tsite_Tsource_reconstruction_quick(fuji_dD, fuji_d18O,method,file_base);
[fuji_T_cond_cor, fuji_T_source_cor, fuji_rs_cor] = Tsite_Tsource_reconstruction_quick(fuji_dD_cor, fuji_d18O_cor,method,file_base);

% [fuji2_T_cond, fuji2_T_source, fuji2_T_surf] = Tsite_Tsource_reconstruction_quick(fuji2_dD, fuji2_d18O,method,file_base);
% [fuji2_T_cond_cor, fuji2_T_source_cor, fuji2_T_surf_cor] = Tsite_Tsource_reconstruction_quick(fuji2_dD_cor, fuji2_d18O_cor,method,file_base);

[SP_T_cond, SP_T_source, SP__rs] = Tsite_Tsource_reconstruction_quick(SP_dD, SP_d18O,method,file_base);
[SP_T_cond_cor, SP_T_source_cor, SP_rs_cor] = Tsite_Tsource_reconstruction_quick(SP_dD_cor, SP_d18O_cor,method,file_base);

[taylor_T_cond, taylor_T_source, taylor_rs] = Tsite_Tsource_reconstruction_quick(taylor_dD, taylor_d18O,method,file_base);
[taylor_T_cond_cor, taylor_T_source_cor, taylor_rs_cor] = Tsite_Tsource_reconstruction_quick(taylor_dD_cor, taylor_d18O_cor,method,file_base);

[talos_T_cond, talos_T_source, talos_rs] = Tsite_Tsource_reconstruction_quick(talos_dD, talos_d18O,method,file_base);
[talos_T_cond_cor, talos_T_source_cor, talos_rs_cor] = Tsite_Tsource_reconstruction_quick(talos_dD_cor, talos_d18O_cor,method,file_base);

%%

figure;
hold on
plot(wais_age,wais_T_cond_cor)
plot(edc_age,edc_T_cond_cor)
plot(edml_age,edml_T_cond_cor)
plot(vostok_age,vostok_T_cond_cor)
plot(siple_age,siple_T_cond_cor)
plot(fuji_age,fuji_T_cond_cor)
plot(SP_age,SP_T_cond_cor,'k')
% plot(taylor_age,taylor_T_cond_cor,'c')
plot(talos_age,talos_T_cond_cor,'c')

figure;
hold on
plot(SP_age,SP_T_cond_cor./0.67,'k')
plot(SP_age,SP_T_cond_cor./0.62,'r')


figure;
hold on
plot(wais_age,wais_T_source_cor)
plot(edc_age,edc_T_source_cor)
plot(edml_age,edml_T_source_cor)
plot(vostok_age,vostok_T_source_cor)
plot(siple_age,siple_T_source_cor)
plot(fuji_age,fuji_T_source_cor)
plot(SP_age,SP_T_source_cor,'k')
% plot(taylor_age,taylor_T_cond_cor,'c')
plot(talos_age,talos_T_source_cor,'c')



figure;
hold on
plot(wais_T_source,wais_T_cond,'.')
plot(edc_T_source,edc_T_cond,'.')
plot(edml_T_source,edml_T_cond,'.')
plot(vostok_T_source,vostok_T_cond,'.')
plot(siple_T_source,siple_T_cond,'.')
plot(fuji_T_source,fuji_T_cond,'.')
plot(SP_T_source,SP_T_cond,'.k')
% plot(taylor_age,taylor_T_cond_cor,'c')
plot(talos_T_source,talos_T_cond,'.c')



figure;
hold on
scatter(wais_T_source,wais_T_cond,10,wais_age)
scatter(edc_T_source,edc_T_cond,10,edc_age)
scatter(edml_T_source,edml_T_cond,10,edml_age)
scatter(vostok_T_source,vostok_T_cond,10,vostok_age)
scatter(siple_T_source,siple_T_cond,10,siple_age)
scatter(fuji_T_source,fuji_T_cond,10,fuji_age)
scatter(SP_T_source,SP_T_cond,10,SP_age)
scatter(talos_T_source,talos_T_cond,10,talos_age)


%%

%% ice core site lats
%from noaa paleoclimate archive
wais_lat=79.467;
edml_lat=75;
edc_lat=75.1;
vostok_lat=78.47;
siple_lat=81.65;
talos_lat=72.8;
fuji_lat=77.32;
SP_lat=90;

%from noaa paleoclimate archive
wais_long=-112.13; %W or E??
edml_long=0.06667;
edc_long=123.4;
vostok_long=106.8;%W
siple_long=208.82;
talos_long=159.1;
fuji_long=38.7;
SP_long=0;

%from NOAA paleocliamte archinve
wais_elev=1766;
edml_elev=2892;
edc_elev=3240;
vostok_elev=3488;
siple_elev=621;
fuji_elev=3810;
talos_elev=2315;
SP_elev=2835;
modern_elev=[wais_elev, edml_elev, edml_elev, vostok_elev, siple_elev, fuji_elev, talos_elev, SP_elev];

%data mostly found on wikipedia
wais_T=-31;
edml_T= -44.6;%from Kohen station
edc_T=-54.5;
vostok_T=-55.18;
siple_T=-25.4;%from NSIDC
fuji_T=-54.3;
talos_T=-41;
SP_T=-49.4;
modern_T=[wais_T, edml_T, edml_T, vostok_T, siple_T, fuji_T, talos_T, SP_T];

%distance to coast from Schoe. 2014, check these!
wais_dist=585;%km
edml_dist=529;%km
edc_dist=870;%km
vostok_dist=1260;%km
siple_dist=470;%km
fuji_dist=1000;%km
talos_dist=250;%km
SP_dist=1400;%km, measured from google earth, about the same from Ross ice shelf or Ronne

Colors=[         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%%
age_cutoff=1000;
figure;
hold on;
plot(wais_T, nanmean(wais_T_cond_cor(wais_age<=age_cutoff)),'o')
plot(edc_T, nanmean(edc_T_cond_cor(edc_age<=age_cutoff)),'o')
plot(edml_T, nanmean(edml_T_cond_cor(edml_age<=age_cutoff)),'o')
plot(SP_T, nanmean(SP_T_cond_cor(SP_age<=age_cutoff)),'ko')
plot(siple_T, nanmean(siple_T_cond_cor(siple_age<=age_cutoff)),'o')
plot(vostok_T, nanmean(vostok_T_cond_cor(vostok_age<=age_cutoff)),'o')
plot(talos_T, nanmean(talos_T_cond_cor(talos_age<=age_cutoff)),'o')
plot([-60:-10],[-60:-10],'k')

%%

wais_T_surf1=(wais_T_cond_cor- -8.1609)./0.6889;%fit MD08 data base assumptions
wais_T_surf2=(wais_T_cond_cor+10.4)./0.67;%fit MD08 data base assumptions, plus cold sources
wais_T_surf3=(wais_T_cond_cor+8.5367)./0.685;%fit MD08 era
wais_T_surf4=(wais_T_cond_cor+15)./0.62;

siple_T_surf1=(siple_T_cond_cor- -8.1609)./0.6889;
siple_T_surf2=(siple_T_cond_cor+10.4)./0.67;
siple_T_surf3=(siple_T_cond_cor+8.5367)./0.685;
siple_T_surf4=(siple_T_cond_cor+15)./0.62;

edc_T_surf1=(edc_T_cond_cor+8.1609)./0.6889;
edc_T_surf2=(edc_T_cond_cor+10.4)./0.67;
edc_T_surf3=(edc_T_cond_cor+8.5367)./0.685;
edc_T_surf4=(edc_T_cond_cor+15)./0.62;

edml_T_surf1=(edml_T_cond_cor+8.1609)./0.6889;
edml_T_surf2=(edml_T_cond_cor+10.4)./0.67;
edml_T_surf3=(edml_T_cond_cor+8.5367)./0.685;
edml_T_surf4=(edml_T_cond_cor+15)./0.62;

talos_T_surf1=(talos_T_cond_cor+8.1609)./0.6889;
talos_T_surf2=(talos_T_cond_cor+10.4)./0.67;
talos_T_surf3=(talos_T_cond_cor+8.5367)./0.685;
talos_T_surf4=(talos_T_cond_cor+15)./0.62;

SP_T_surf1=(SP_T_cond_cor- -8.1609)./0.6889;
SP_T_surf2=(SP_T_cond_cor+10.4)./0.67;
SP_T_surf3=(SP_T_cond_cor+8.5367)./0.685;
SP_T_surf4=(SP_T_cond_cor+15)./0.62;

vostok_T_surf1=(vostok_T_cond_cor- -8.1609)./0.6889;
vostok_T_surf2=(vostok_T_cond_cor+10.4)./0.67;
vostok_T_surf3=(vostok_T_cond_cor+8.5367)./0.685;
vostok_T_surf4=(vostok_T_cond_cor+15)./0.62;

fuji_T_surf1=(fuji_T_cond_cor- -8.1609)./0.6889;
fuji_T_surf2=(fuji_T_cond_cor+10.4)./0.67;
fuji_T_surf3=(fuji_T_cond_cor+8.5367)./0.685;
fuji_T_surf4=(fuji_T_cond_cor+15)./0.62;


figure;
hold on
plot(wais_age,wais_T_surf1,'k')
plot(wais_age,wais_T_surf2,'r')
plot(wais_age,wais_T_surf3,'b')
plot(wais_age,wais_T_surf4,'c')

plot(edc_age,edc_T_surf1,'k')
plot(edc_age,edc_T_surf2,'r')
plot(edc_age,edc_T_surf3,'b')
plot(edc_age,edc_T_surf4,'c')

plot(SP_age,SP_T_surf1,'k')
plot(SP_age,SP_T_surf2,'r')
plot(SP_age,SP_T_surf3,'b')
plot(SP_age,SP_T_surf4,'c')

plot(vostok_age,vostok_T_surf1,'k')
plot(vostok_age,vostok_T_surf2,'r')
plot(vostok_age,vostok_T_surf3,'b')
plot(vostok_age,vostok_T_surf4,'c')



 

cutoff=50;
figure(10);
hold on
plot([-60 -20],[-60 -20],'k')
scatter(wais_T,nanmean(wais_T_surf1(wais_age<cutoff)),'filled','k')
scatter(SP_T,nanmean(SP_T_surf1(SP_age<cutoff)),'filled','k')
scatter(edc_T,nanmean(edc_T_surf1(edc_age<cutoff)),'filled','k')
scatter(edml_T,nanmean(edml_T_surf1(edml_age<cutoff)),'filled','k')
scatter(vostok_T,nanmean(vostok_T_surf1(vostok_age<cutoff)),'filled','k')
scatter(siple_T,nanmean(siple_T_surf1(siple_age<cutoff)),'filled','k')
scatter(talos_T,nanmean(talos_T_surf1(talos_age<cutoff)),'filled','k')

scatter(wais_T,nanmean(wais_T_surf2(wais_age<cutoff)),'filled','r')
scatter(SP_T,nanmean(SP_T_surf2(SP_age<cutoff)),'filled','r')
scatter(edc_T,nanmean(edc_T_surf2(edc_age<cutoff)),'filled','r')
scatter(edml_T,nanmean(edml_T_surf2(edml_age<cutoff)),'filled','r')
scatter(vostok_T,nanmean(vostok_T_surf2(vostok_age<cutoff)),'filled','r')
scatter(siple_T,nanmean(siple_T_surf2(siple_age<cutoff)),'filled','r')
scatter(talos_T,nanmean(talos_T_surf2(talos_age<cutoff)),'filled','r')

scatter(wais_T,nanmean(wais_T_surf3(wais_age<cutoff)),'filled','b')
scatter(SP_T,nanmean(SP_T_surf3(SP_age<cutoff)),'filled','b')
scatter(edc_T,nanmean(edc_T_surf3(edc_age<cutoff)),'filled','b')
scatter(edml_T,nanmean(edml_T_surf3(edml_age<cutoff)),'filled','b')
scatter(vostok_T,nanmean(vostok_T_surf3(vostok_age<cutoff)),'filled','b')
scatter(siple_T,nanmean(siple_T_surf3(siple_age<cutoff)),'filled','b')
scatter(talos_T,nanmean(talos_T_surf3(talos_age<cutoff)),'filled','b')

scatter(wais_T,nanmean(wais_T_surf4(wais_age<cutoff)),'filled','c')
scatter(SP_T,nanmean(SP_T_surf4(SP_age<cutoff)),'filled','c')
scatter(edc_T,nanmean(edc_T_surf4(edc_age<cutoff)),'filled','c')
scatter(edml_T,nanmean(edml_T_surf4(edml_age<cutoff)),'filled','c')
scatter(vostok_T,nanmean(vostok_T_surf4(vostok_age<cutoff)),'filled','c')
scatter(siple_T,nanmean(siple_T_surf4(siple_age<cutoff)),'filled','c')
scatter(talos_T,nanmean(talos_T_surf4(talos_age<cutoff)),'filled','c')


XXX=[wais_T SP_T edc_T edml_T vostok_T siple_T talos_T];
YYY = [nanmean(wais_T_surf1(wais_age<cutoff)) nanmean(SP_T_surf1(SP_age<cutoff)) nanmean(edc_T_surf1(edc_age<cutoff)) nanmean(edml_T_surf1(edml_age<cutoff)) ...
    nanmean(vostok_T_surf1(vostok_age<cutoff)) nanmean(siple_T_surf1(siple_age<cutoff)) nanmean(talos_T_surf1(talos_age<cutoff))];
idx=find(~isnan(XXX) & ~isnan(YYY));
P_surf1 = polyfit(XXX(idx),YYY(idx),1) 


XXX=[wais_T SP_T edc_T edml_T vostok_T siple_T talos_T];
YYY = [nanmean(wais_T_surf2(wais_age<cutoff)) nanmean(SP_T_surf2(SP_age<cutoff)) nanmean(edc_T_surf2(edc_age<cutoff)) nanmean(edml_T_surf2(edml_age<cutoff)) ...
    nanmean(vostok_T_surf2(vostok_age<cutoff)) nanmean(siple_T_surf2(siple_age<cutoff)) nanmean(talos_T_surf2(talos_age<cutoff))];
idx=find(~isnan(XXX) & ~isnan(YYY));
P_surf2 = polyfit(XXX(idx),YYY(idx),1) 

XXX=[wais_T SP_T edc_T edml_T vostok_T siple_T talos_T];
YYY = [nanmean(wais_T_surf3(wais_age<cutoff)) nanmean(SP_T_surf3(SP_age<cutoff)) nanmean(edc_T_surf3(edc_age<cutoff)) nanmean(edml_T_surf3(edml_age<cutoff)) ...
    nanmean(vostok_T_surf3(vostok_age<cutoff)) nanmean(siple_T_surf3(siple_age<cutoff)) nanmean(talos_T_surf3(talos_age<cutoff))];
idx=find(~isnan(XXX) & ~isnan(YYY));
P_surf3 = polyfit(XXX(idx),YYY(idx),1) 



XXX=[wais_T SP_T edc_T edml_T vostok_T siple_T talos_T];
YYY = [nanmean(wais_T_cond_cor(wais_age<cutoff)) nanmean(SP_T_cond_cor(SP_age<cutoff)) nanmean(edc_T_cond_cor(edc_age<cutoff)) nanmean(edml_T_cond_cor(edml_age<cutoff)) ...
    nanmean(vostok_T_cond_cor(vostok_age<cutoff)) nanmean(siple_T_cond_cor(siple_age<cutoff)) nanmean(talos_T_cond_cor(talos_age<cutoff))];
idx=find(~isnan(XXX) & ~isnan(YYY));
P_cond = polyfit(XXX(idx),YYY(idx),1) 


XXX=[wais_T SP_T edc_T edml_T vostok_T siple_T talos_T];
cutoff=50;
YYY_50 = [nanmean(wais_T_cond(wais_age<cutoff)) nanmean(SP_T_cond(SP_age<cutoff)) nanmean(edc_T_cond(edc_age<cutoff)) nanmean(edml_T_cond(edml_age<cutoff)) ...
    nanmean(vostok_T_cond(vostok_age<cutoff)) nanmean(siple_T_cond(siple_age<cutoff)) nanmean(talos_T_cond(talos_age<cutoff))];
idx_50=find(~isnan(XXX) & ~isnan(YYY_50));
cutoff=50;
YYY_100 = [nanmean(wais_T_cond(wais_age<cutoff)) nanmean(SP_T_cond(SP_age<cutoff)) nanmean(edc_T_cond(edc_age<cutoff)) nanmean(edml_T_cond(edml_age<cutoff)) ...
    nanmean(vostok_T_cond(vostok_age<cutoff)) nanmean(siple_T_cond(siple_age<cutoff)) nanmean(talos_T_cond(talos_age<cutoff))];
idx_100=find(~isnan(XXX) & ~isnan(YYY_100));
ice_core_Ts=XXX;
ice_core_Tc_50=YYY_50;
ice_core_Tc_100=YYY_100;
save('ice_core_Ts_Tc.mat','ice_core_Ts','ice_core_Tc_50','ice_core_Tc_100','idx_50','idx_100');

figure;
hold on
plot(wais_age,wais_T_surf1)
plot(siple_age,siple_T_surf1)
plot(edc_age,edc_T_surf1)
plot(edml_age,edml_T_surf1)
plot(SP_age,SP_T_surf1)
plot(vostok_age,vostok_T_surf1)
plot(fuji_age,fuji_T_surf1)
plot(talos_age,talos_T_surf1)


% LGM_age=[19000 21000];
LGM_age=[20000 30000];

% EH_age=[9000 11000];
EH_age=[8000 10000];

LH_age=[0 10000];

figure;
hold on
plot(wais_T.*ones(size(wais_T_surf1)),wais_T_surf1)

figure;
hold on
scatter(wais_T,nanmean(wais_T_surf1(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
scatter(siple_T,nanmean(siple_T_surf1(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
scatter(SP_T,nanmean(SP_T_surf1(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled');
scatter(edc_T,nanmean(edc_T_surf1(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
scatter(edml_T,nanmean(edml_T_surf1(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
scatter(fuji_T,nanmean(fuji_T_surf1(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
scatter(vostok_T,nanmean(vostok_T_surf1(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
scatter(talos_T,nanmean(talos_T_surf1(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');



figure;
hold on
scatter(wais_T,nanmean(wais_T_surf1(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled','k');
scatter(siple_T,nanmean(siple_T_surf1(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled','k');
scatter(SP_T,nanmean(SP_T_surf1(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled','k');
scatter(edc_T,nanmean(edc_T_surf1(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled','k');
scatter(edml_T,nanmean(edml_T_surf1(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled','k');
scatter(fuji_T,nanmean(fuji_T_surf1(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled','k');
scatter(vostok_T,nanmean(vostok_T_surf1(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled','k');
scatter(talos_T,nanmean(talos_T_surf1(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled','k');


scatter(wais_T,nanmean(wais_T_surf2(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_surf2(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled','r');
scatter(siple_T,nanmean(siple_T_surf2(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_surf2(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled','r');
scatter(SP_T,nanmean(SP_T_surf2(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_surf2(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled','r');
scatter(edc_T,nanmean(edc_T_surf2(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_surf2(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled','r');
scatter(edml_T,nanmean(edml_T_surf2(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_surf2(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled','r');
scatter(fuji_T,nanmean(fuji_T_surf2(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_surf2(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled','r');
scatter(vostok_T,nanmean(vostok_T_surf2(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_surf2(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled','r');
scatter(talos_T,nanmean(talos_T_surf2(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_surf2(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled','r');

scatter(wais_T,nanmean(wais_T_surf3(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_surf3(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled','b');
scatter(siple_T,nanmean(siple_T_surf3(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_surf3(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled','b');
scatter(SP_T,nanmean(SP_T_surf3(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_surf3(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled','b');
scatter(edc_T,nanmean(edc_T_surf3(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_surf3(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled','b');
scatter(edml_T,nanmean(edml_T_surf3(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_surf3(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled','b');
scatter(fuji_T,nanmean(fuji_T_surf3(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_surf3(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled','b');
scatter(vostok_T,nanmean(vostok_T_surf3(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_surf3(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled','b');
scatter(talos_T,nanmean(talos_T_surf3(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_surf3(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled','b');


figure;
hold on
scatter(wais_T,nanmean(wais_T_surf1(wais_age>=LH_age(1) & wais_age<=LH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
scatter(siple_T,nanmean(siple_T_surf1(siple_age>=LH_age(1) & siple_age<=LH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
scatter(SP_T,nanmean(SP_T_surf1(SP_age>=LH_age(1) & SP_age<=LH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled','k');
scatter(edc_T,nanmean(edc_T_surf1(edc_age>=LH_age(1) & edc_age<=LH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
scatter(edml_T,nanmean(edml_T_surf1(edml_age>=LH_age(1) & edml_age<=LH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
scatter(fuji_T,nanmean(fuji_T_surf1(fuji_age>=LH_age(1) & fuji_age<=LH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
scatter(vostok_T,nanmean(vostok_T_surf1(vostok_age>=LH_age(1) & vostok_age<=LH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
scatter(talos_T,nanmean(talos_T_surf1(talos_age>=LH_age(1) & talos_age<=LH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');


%%
load 'msd_lats_ice_core_sites_3.mat'


figure;
hold on
scatter(-wais_lat,nanmean(wais_T_surf1(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
scatter(-siple_lat,nanmean(siple_T_surf1(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
scatter(-SP_lat,nanmean(SP_T_surf1(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled');
scatter(-edc_lat,nanmean(edc_T_surf1(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
scatter(-edml_lat,nanmean(edml_T_surf1(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
scatter(-fuji_lat,nanmean(fuji_T_surf1(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
scatter(-vostok_lat,nanmean(vostok_T_surf1(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
scatter(-talos_lat,nanmean(talos_T_surf1(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');



scatter(lati_wais,nanmean(wais_T_source_cor(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_source_cor(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
scatter(lati_siple,nanmean(siple_T_source_cor(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_source_cor(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
scatter(lati_SP,nanmean(SP_T_source_cor(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_source_cor(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled');
scatter(lati_edc,nanmean(edc_T_source_cor(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_source_cor(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
scatter(lati_edml,nanmean(edml_T_source_cor(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_source_cor(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
scatter(lati_fuji,nanmean(fuji_T_source_cor(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_source_cor(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
scatter(lati_vostok,nanmean(vostok_T_source_cor(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_source_cor(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
scatter(lati_talos,nanmean(talos_T_source_cor(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_source_cor(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');


%% Margo stuff
%MARGO defines LGM at 19-23ka
LGM_age=[19000 23000];
%for there anomalies MARGO uses modern. soo...?
LH_age=[0 4000];
LH_age2=[0 1000];


ncdisp('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc')
lon=ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','lon');
lat=ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','lat');

[coast]=load('coast');
coast.long(coast.long<=0)=coast.long(coast.long<=0)+360;


%'WOA (1998) 10-m temperature (annual)'
woaann = ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','woaann');
%'MARGO (2009) reconstructed LGM SST (annual)'
recann = ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','recann');
%'MARGO (2009) LGM SST anomaly (annual)'
anoann = ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','anoann');
%'MARGO (2009) reconstructed LGM SST (NH summer/SH winter)'
 recjas = ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','recjas'); 
 % 'MARGO (2009) reconstructed LGM SST (NH winter/SH summer)'
 recjfm = ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','recjfm');  
 %  'MARGO (2009) LGM SST anomaly (NH summer/SH winter)'
 anojas = ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','anojas');
 %  'MARGO (2009) LGM SST anomaly (NH winter/SH summer)'
anojfm = ncread('/Users/bradley/work/2020/LGM_SSTs/MARGO/margo_all_lgm_sst_with_errors_weighted_mean.nc','anojfm'); 


%zonal means
anoann_zonal=nanmean(anoann);
anojfm_zonal=nanmean(anojfm);
anojas_zonal=nanmean(anojas);

fig('units','inches','width',8,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(lat,-anoann,'ko');
plot(lat,-anoann_zonal,'Linewidth',2);
xlim([-90 10])
% ylim([-10 5])
title('ANN')


fig('units','inches','width',8,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(lat,recjfm,'ko');
xlim([-90 10])


%%
fig(500,'units','inches','width',8,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
xlim([-90 0])

plot(lat,-anoann,'ko');

% 
% scatter(-wais_lat,nanmean(wais_T_surf1(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
% scatter(-siple_lat,nanmean(siple_T_surf1(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
% scatter(-SP_lat,nanmean(SP_T_surf1(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled');
% scatter(-edc_lat,nanmean(edc_T_surf1(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
% scatter(-edml_lat,nanmean(edml_T_surf1(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
% scatter(-fuji_lat,nanmean(fuji_T_surf1(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
% scatter(-vostok_lat,nanmean(vostok_T_surf1(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
% scatter(-talos_lat,nanmean(talos_T_surf1(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');
% 
% scatter(lati_wais,nanmean(wais_T_source_cor(wais_age>=EH_age(1) & wais_age<=EH_age(2)))-nanmean(wais_T_source_cor(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
% scatter(lati_siple,nanmean(siple_T_source_cor(siple_age>=EH_age(1) & siple_age<=EH_age(2)))-nanmean(siple_T_source_cor(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
% scatter(lati_SP,nanmean(SP_T_source_cor(SP_age>=EH_age(1) & SP_age<=EH_age(2)))-nanmean(SP_T_source_cor(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled');
% scatter(lati_edc,nanmean(edc_T_source_cor(edc_age>=EH_age(1) & edc_age<=EH_age(2)))-nanmean(edc_T_source_cor(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
% scatter(lati_edml,nanmean(edml_T_source_cor(edml_age>=EH_age(1) & edml_age<=EH_age(2)))-nanmean(edml_T_source_cor(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
% scatter(lati_fuji,nanmean(fuji_T_source_cor(fuji_age>=EH_age(1) & fuji_age<=EH_age(2)))-nanmean(fuji_T_source_cor(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
% scatter(lati_vostok,nanmean(vostok_T_source_cor(vostok_age>=EH_age(1) & vostok_age<=EH_age(2)))-nanmean(vostok_T_source_cor(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
% scatter(lati_talos,nanmean(talos_T_source_cor(talos_age>=EH_age(1) & talos_age<=EH_age(2)))-nanmean(talos_T_source_cor(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');

scatter(-wais_lat,nanmean(wais_T_surf1(wais_age>=LH_age(1) & wais_age<=LH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
scatter(-siple_lat,nanmean(siple_T_surf1(siple_age>=LH_age(1) & siple_age<=LH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
scatter(-SP_lat,nanmean(SP_T_surf1(SP_age>=LH_age(1) & SP_age<=LH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled');
scatter(-edc_lat,nanmean(edc_T_surf1(edc_age>=LH_age(1) & edc_age<=LH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
scatter(-edml_lat,nanmean(edml_T_surf1(edml_age>=LH_age(1) & edml_age<=LH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
scatter(-fuji_lat,nanmean(fuji_T_surf1(fuji_age>=LH_age(1) & fuji_age<=LH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
scatter(-vostok_lat,nanmean(vostok_T_surf1(vostok_age>=LH_age(1) & vostok_age<=LH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
scatter(-talos_lat,nanmean(talos_T_surf1(talos_age>=LH_age(1) & talos_age<=LH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');

scatter(lati_wais,nanmean(wais_T_source_cor(wais_age>=LH_age(1) & wais_age<=LH_age(2)))-nanmean(wais_T_source_cor(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled');
scatter(lati_siple,nanmean(siple_T_source_cor(siple_age>=LH_age(1) & siple_age<=LH_age(2)))-nanmean(siple_T_source_cor(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled');
scatter(lati_SP,nanmean(SP_T_source_cor(SP_age>=LH_age(1) & SP_age<=LH_age(2)))-nanmean(SP_T_source_cor(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled');
scatter(lati_edc,nanmean(edc_T_source_cor(edc_age>=LH_age(1) & edc_age<=LH_age(2)))-nanmean(edc_T_source_cor(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled');
scatter(lati_edml,nanmean(edml_T_source_cor(edml_age>=LH_age(1) & edml_age<=LH_age(2)))-nanmean(edml_T_source_cor(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled');
scatter(lati_fuji,nanmean(fuji_T_source_cor(fuji_age>=LH_age(1) & fuji_age<=LH_age(2)))-nanmean(fuji_T_source_cor(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled');
scatter(lati_vostok,nanmean(vostok_T_source_cor(vostok_age>=LH_age(1) & vostok_age<=LH_age(2)))-nanmean(vostok_T_source_cor(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled');
scatter(lati_talos,nanmean(talos_T_source_cor(talos_age>=LH_age(1) & talos_age<=LH_age(2)))-nanmean(talos_T_source_cor(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled');


%%

SZ=70;
figure;
hold on
xlim([-90 0])

plot(lat,recann ,'ko');
scatter(lati_wais,nanmean(wais_T_source_cor(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),SZ,'filled');
scatter(lati_siple,nanmean(siple_T_source_cor(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),SZ,'filled');
scatter(lati_SP,nanmean(SP_T_source_cor(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),SZ,'filled');
scatter(lati_edc,nanmean(edc_T_source_cor(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),SZ,'filled');
scatter(lati_edml,nanmean(edml_T_source_cor(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),SZ,'filled');
scatter(lati_fuji,nanmean(fuji_T_source_cor(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),SZ,'filled');
scatter(lati_vostok,nanmean(vostok_T_source_cor(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),SZ,'filled');
scatter(lati_talos,nanmean(talos_T_source_cor(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),SZ,'filled');



figure;
hold on
xlim([-90 0])

plot(lat,woaann ,'ko');
scatter(lati_wais,nanmean(wais_T_source_cor(wais_age>=LH_age(1) & wais_age<=LH_age(2))),'filled');
scatter(lati_siple,nanmean(siple_T_source_cor(siple_age>=LH_age(1) & siple_age<=LH_age(2))),'filled');
scatter(lati_SP,nanmean(SP_T_source_cor(SP_age>=LH_age(1) & SP_age<=LH_age(2))),'filled');
scatter(lati_edc,nanmean(edc_T_source_cor(edc_age>=LH_age(1) & edc_age<=LH_age(2))),'filled');
scatter(lati_edml,nanmean(edml_T_source_cor(edml_age>=LH_age(1) & edml_age<=LH_age(2))),'filled');
scatter(lati_fuji,nanmean(fuji_T_source_cor(fuji_age>=LH_age(1) & fuji_age<=LH_age(2))),'filled');
scatter(lati_vostok,nanmean(vostok_T_source_cor(vostok_age>=LH_age(1) & vostok_age<=LH_age(2))),'filled');
scatter(lati_talos,nanmean(talos_T_source_cor(talos_age>=LH_age(1) & talos_age<=LH_age(2))),'filled');


%%

[rh0, delta_rh0, sst0, delta_sst0, RHn0, deltarhn0] = T_RH_RHn_2020(0:25, 1, 'ncep', 'spline');

T_source_LH=[nanmean(wais_T_source_cor(wais_age>=LH_age(1) & wais_age<=LH_age(2))),...
    nanmean(siple_T_source_cor(siple_age>=LH_age(1) & siple_age<=LH_age(2))),...
    nanmean(SP_T_source_cor(SP_age>=LH_age(1) & SP_age<=LH_age(2))),...
    nanmean(edc_T_source_cor(edc_age>=LH_age(1) & edc_age<=LH_age(2))),...
    nanmean(edml_T_source_cor(edml_age>=LH_age(1) & edml_age<=LH_age(2))),...
    nanmean(fuji_T_source_cor(fuji_age>=LH_age(1) & fuji_age<=LH_age(2))),...
    nanmean(vostok_T_source_cor(vostok_age>=LH_age(1) & vostok_age<=LH_age(2))),...
    nanmean(talos_T_source_cor(talos_age>=LH_age(1) & talos_age<=LH_age(2))),...
    ];

T_source_LH2=[nanmean(wais_T_source_cor(wais_age>=LH_age2(1) & wais_age<=LH_age2(2))),...
    nanmean(siple_T_source_cor(siple_age>=LH_age2(1) & siple_age<=LH_age2(2))),...
    nanmean(SP_T_source_cor(SP_age>=LH_age2(1) & SP_age<=LH_age2(2))),...
    nanmean(edc_T_source_cor(edc_age>=LH_age2(1) & edc_age<=LH_age2(2))),...
    nanmean(edml_T_source_cor(edml_age>=LH_age2(1) & edml_age<=LH_age2(2))),...
    nanmean(fuji_T_source_cor(fuji_age>=LH_age2(1) & fuji_age<=LH_age2(2))),...
    nanmean(vostok_T_source_cor(vostok_age>=LH_age2(1) & vostok_age<=LH_age2(2))),...
    nanmean(talos_T_source_cor(talos_age>=LH_age2(1) & talos_age<=LH_age2(2))),...
    ];

SST_source_LH=interp1(0:25,sst0,T_source_LH);
SST_source_LH2=interp1(0:25,sst0,T_source_LH2);

T_source_LGM=[nanmean(wais_T_source_cor(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),...
    nanmean(siple_T_source_cor(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),...
    nanmean(SP_T_source_cor(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),...
    nanmean(edc_T_source_cor(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),...
    nanmean(edml_T_source_cor(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),...
    nanmean(fuji_T_source_cor(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),...
    nanmean(vostok_T_source_cor(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),...
    nanmean(talos_T_source_cor(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),...
    ];

SST_source_LGM=interp1(0:25,sst0,T_source_LGM);

lati_all=[lati_wais, lati_siple, lati_SP, lati_edc, lati_edml, lati_fuji, lati_vostok, lati_talos];

figure;
hold on
xlim([-90 0])
plot(lat,woaann ,'ko');
scatter(lati_all,SST_source_LH)



figure;
hold on
xlim([-90 0])
plot(lat,recann ,'ko');
scatter(lati_all,SST_source_LGM)

[lats,lons] = meshgrid(lat,lon);

fig('units','inches','width',8,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
xlim([-90 0])

scatter(lats(:),-anoann(:));
scatter(-wais_lat,nanmean(wais_T_surf1(wais_age>=LH_age(1) & wais_age<=LH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'b','filled');
scatter(-siple_lat,nanmean(siple_T_surf1(siple_age>=LH_age(1) & siple_age<=LH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'b','filled');
scatter(-SP_lat,nanmean(SP_T_surf1(SP_age>=LH_age(1) & SP_age<=LH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'b','filled');
scatter(-edc_lat,nanmean(edc_T_surf1(edc_age>=LH_age(1) & edc_age<=LH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'b','filled');
scatter(-edml_lat,nanmean(edml_T_surf1(edml_age>=LH_age(1) & edml_age<=LH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'b','filled');
scatter(-fuji_lat,nanmean(fuji_T_surf1(fuji_age>=LH_age(1) & fuji_age<=LH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'b','filled');
scatter(-vostok_lat,nanmean(vostok_T_surf1(vostok_age>=LH_age(1) & vostok_age<=LH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'b','filled');
scatter(-talos_lat,nanmean(talos_T_surf1(talos_age>=LH_age(1) & talos_age<=LH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'b','filled');


scatter(-wais_lat,nanmean(wais_T_surf1(wais_age>=LH_age2(1) & wais_age<=LH_age2(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'b');
scatter(-siple_lat,nanmean(siple_T_surf1(siple_age>=LH_age2(1) & siple_age<=LH_age2(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'b');
scatter(-SP_lat,nanmean(SP_T_surf1(SP_age>=LH_age2(1) & SP_age<=LH_age2(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'b');
scatter(-edc_lat,nanmean(edc_T_surf1(edc_age>=LH_age2(1) & edc_age<=LH_age2(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'b');
scatter(-edml_lat,nanmean(edml_T_surf1(edml_age>=LH_age2(1) & edml_age<=LH_age2(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'b');
scatter(-fuji_lat,nanmean(fuji_T_surf1(fuji_age>=LH_age2(1) & fuji_age<=LH_age2(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'b');
scatter(-vostok_lat,nanmean(vostok_T_surf1(vostok_age>=LH_age2(1) & vostok_age<=LH_age2(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'b');
scatter(-talos_lat,nanmean(talos_T_surf1(talos_age>=LH_age2(1) & talos_age<=LH_age2(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'b');


scatter(lati_all,SST_source_LH-SST_source_LGM,'filled')
scatter(lati_all,SST_source_LH2-SST_source_LGM)



%%
longi_all=[wais_long, siple_long, SP_long, edc_long, edml_long, fuji_long, vostok_long, talos_long];

load coastlines


fig('units','inches','width',8,'height',8,'font','Helvetica','fontsize',16,'border','on');
hold on
% worldmap('World')
worldmap([-90 -2],[0 360])

plotm(coastlat,coastlon,'k')

scatterm(lats(:),lons(:),20,-anoann(:),'filled');

scatterm(lati_all,longi_all,100,SST_source_LH-SST_source_LGM,'filled');


scatterm(-wais_lat,wais_long,100,nanmean(wais_T_surf1(wais_age>=LH_age(1) & wais_age<=LH_age(2)))-nanmean(wais_T_surf1(wais_age>=LGM_age(1) & wais_age<=LGM_age(2))),'filled')
scatterm(-siple_lat,siple_long,100,nanmean(siple_T_surf1(siple_age>=LH_age(1) & siple_age<=LH_age(2)))-nanmean(siple_T_surf1(siple_age>=LGM_age(1) & siple_age<=LGM_age(2))),'filled')
scatterm(-SP_lat,SP_long,100,nanmean(SP_T_surf1(SP_age>=LH_age(1) & SP_age<=LH_age(2)))-nanmean(SP_T_surf1(SP_age>=LGM_age(1) & SP_age<=LGM_age(2))),'filled')
scatterm(-edc_lat,edc_long,100,nanmean(edc_T_surf1(edc_age>=LH_age(1) & edc_age<=LH_age(2)))-nanmean(edc_T_surf1(edc_age>=LGM_age(1) & edc_age<=LGM_age(2))),'filled')
scatterm(-edml_lat,edml_long,100,nanmean(edml_T_surf1(edml_age>=LH_age(1) & edml_age<=LH_age(2)))-nanmean(edml_T_surf1(edml_age>=LGM_age(1) & edml_age<=LGM_age(2))),'filled')
scatterm(-vostok_lat,vostok_long,100,nanmean(vostok_T_surf1(vostok_age>=LH_age(1) & vostok_age<=LH_age(2)))-nanmean(vostok_T_surf1(vostok_age>=LGM_age(1) & vostok_age<=LGM_age(2))),'filled')
scatterm(-fuji_lat,fuji_long,100,nanmean(fuji_T_surf1(fuji_age>=LH_age(1) & fuji_age<=LH_age(2)))-nanmean(fuji_T_surf1(fuji_age>=LGM_age(1) & fuji_age<=LGM_age(2))),'filled')
scatterm(-talos_lat,talos_long,100,nanmean(talos_T_surf1(talos_age>=LH_age(1) & talos_age<=LH_age(2)))-nanmean(talos_T_surf1(talos_age>=LGM_age(1) & talos_age<=LGM_age(2))),'filled')

%% 

% d18Osw_i=-0.3;
% [wais_d18O_ln_cor,wais_dD_ln_cor,wais_dln_cor, wais_d18O_cor, wais_dD_cor]= seawater_cor_ln(wais_d18O,wais_dD,wais_age',d18Osw_i);
% 
% 
% file_base_hi = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_0_28_SITE_m70_0.mat';
% file_base_lo = './2020/SWIM_results_202006_ncep_annual_adiabat_local_initsat_a1_b00525_c0_adiffevap_009_SOURCE_m10_0_SITE_m70_0.mat';
% 
% [wais_T_cond, wais_T_source, wais_rs] = Tsite_Tsource_reconstruction_quick(wais_dD, wais_d18O,method,file_base);
% [wais_T_cond_cor, wais_T_source_cor, wais_rs_cor] = Tsite_Tsource_reconstruction_quick(wais_dD_cor, wais_d18O_cor,method,file_base);
% wdc_T_surf_base=Ts_to_Tc_2020(wdc_T_cond_base,1);
