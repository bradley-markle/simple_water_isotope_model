%%


%% 
Colors=[         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%% load ship vapor data

load('Ship_ACTIV.csv');
load('Ship_JARE55.csv');
load('Ship_STRASSE.csv');
load('Ship_RARAAVIS.csv');
load('Ship_Bermuda.csv');

Ship_ACTIV(Ship_ACTIV==-9999)=NaN;
Ship_JARE55(Ship_JARE55==-9999)=NaN;
Ship_STRASSE(Ship_STRASSE==-9999)=NaN;
Ship_RARAAVIS(Ship_RARAAVIS==-9999)=NaN;
Ship_Bermuda(Ship_Bermuda==-9999)=NaN;

Ship_ACTIV(Ship_ACTIV==-99.99)=NaN;
Ship_JARE55(Ship_JARE55==-99.99)=NaN;
Ship_STRASSE(Ship_STRASSE==-99.99)=NaN;
Ship_RARAAVIS(Ship_RARAAVIS==-99.99)=NaN;
Ship_Bermuda(Ship_Bermuda==-99.99)=NaN;

latitude_ACTIV=Ship_ACTIV(:,2);
longitude_ACTIV=Ship_ACTIV(:,3);
d18O_ACTIV=Ship_ACTIV(:,5);
dD_ACTIV=Ship_ACTIV(:,7);
Ta_ACTIV=Ship_ACTIV(:,9);
RH_ACTIV=Ship_ACTIV(:,10);

latitude_JARE55=Ship_JARE55(:,2);
longitude_JARE55=Ship_JARE55(:,3);
d18O_JARE55=Ship_JARE55(:,5);
dD_JARE55=Ship_JARE55(:,7);
Ta_JARE55=Ship_JARE55(:,9);
RH_JARE55=Ship_JARE55(:,10);

latitude_STRASSE=Ship_STRASSE(:,2);
longitude_STRASSE=Ship_STRASSE(:,3);
d18O_STRASSE=Ship_STRASSE(:,5);
dD_STRASSE=Ship_STRASSE(:,7);
Ta_STRASSE=Ship_STRASSE(:,9);
RH_STRASSE=Ship_STRASSE(:,10);

latitude_RARAAVIS=Ship_RARAAVIS(:,2);
longitude_RARAAVIS=Ship_RARAAVIS(:,3);
d18O_RARAAVIS=Ship_RARAAVIS(:,5);
dD_RARAAVIS=Ship_RARAAVIS(:,7);
Ta_RARAAVIS=Ship_RARAAVIS(:,9);
RH_RARAAVIS=Ship_RARAAVIS(:,10);

latitude_Bermuda=Ship_Bermuda(:,2);
longitude_Bermuda=Ship_Bermuda(:,3);
d18O_Bermuda=Ship_Bermuda(:,5);
dD_Bermuda=Ship_Bermuda(:,7);
Ta_Bermuda=Ship_Bermuda(:,9);
RH_Bermuda=Ship_Bermuda(:,10);


fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(latitude_JARE55,Ta_JARE55,'.')
plot(latitude_STRASSE,Ta_STRASSE,'.')
plot(latitude_RARAAVIS,Ta_RARAAVIS,'.')
plot(latitude_ACTIV,Ta_ACTIV,'.')
plot(latitude_Bermuda,Ta_Bermuda,'.')


fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(latitude_JARE55,RH_JARE55,'.')
plot(latitude_STRASSE,RH_STRASSE,'.')
plot(latitude_RARAAVIS,RH_RARAAVIS,'.')
plot(latitude_ACTIV,RH_ACTIV,'.')
plot(latitude_Bermuda,RH_Bermuda,'.')



fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(latitude_JARE55,d18O_JARE55,'.')
plot(latitude_STRASSE,d18O_STRASSE,'.')
plot(latitude_RARAAVIS,d18O_RARAAVIS,'.')
plot(latitude_ACTIV,d18O_ACTIV,'.')
plot(latitude_Bermuda,d18O_Bermuda,'.')

fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(Ta_JARE55,d18O_JARE55,'.')
plot(Ta_ACTIV,d18O_ACTIV,'.')
plot(Ta_STRASSE,d18O_STRASSE,'.')
plot(Ta_RARAAVIS,d18O_RARAAVIS,'.')
plot(Ta_Bermuda,d18O_Bermuda,'.')

fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(Ta_JARE55,dD_JARE55,'.')
plot(Ta_ACTIV,dD_ACTIV,'.')
plot(Ta_STRASSE,dD_STRASSE,'.')
plot(Ta_RARAAVIS,dD_RARAAVIS,'.')
plot(Ta_Bermuda,dD_Bermuda,'.')




fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(RH_JARE55,d18O_JARE55,'.')
plot(RH_ACTIV,d18O_ACTIV,'.')
plot(RH_STRASSE,d18O_STRASSE,'.')
plot(RH_RARAAVIS,d18O_RARAAVIS,'.')
plot(RH_Bermuda,d18O_Bermuda,'.')


%% calculate ssts
cd ../
[~, ~, sst_ncep_ACTIV, ~, ~, ~] = T_RH_RHn_NCEP(Ta_ACTIV, 0);
[~, ~, sst_era_ACTIV, ~, ~, ~, ~, ~] = T_RH_RHn_ERA(Ta_ACTIV, 0, 'annual');

[~, ~, sst_ncep_JARE55, ~, ~, ~] = T_RH_RHn_NCEP(Ta_JARE55, 1);
[~, ~, sst_era_JARE55, ~, ~, ~, ~, ~] = T_RH_RHn_ERA(Ta_JARE55, 1, 'annual');
[~, ~, sst_era_alt_JARE55, ~, ~, ~, ~, ~] = T_RH_RHn_ERA(Ta_JARE55, 1, 'DJF');

[~, ~, sst_ncep_STRASSE, ~, ~, ~] = T_RH_RHn_NCEP(Ta_STRASSE, 0);
[~, ~, sst_era_STRASSE, ~, ~, ~, ~, ~] = T_RH_RHn_ERA(Ta_STRASSE, 0, 'annual');

[~, ~, sst_ncep_RARAAVIS, ~, ~, ~] = T_RH_RHn_NCEP(Ta_RARAAVIS, 0);
[~, ~, sst_era_RARAAVIS, ~, ~, ~, ~, ~] = T_RH_RHn_ERA(Ta_RARAAVIS, 0, 'annual');

[~, ~, sst_ncep_Bermuda, ~, ~, ~] = T_RH_RHn_NCEP(Ta_Bermuda, 0);
[~, ~, sst_era_Bermuda, ~, ~, ~, ~, ~] = T_RH_RHn_ERA(Ta_Bermuda, 0, 'annual');


fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(sst_ncep_JARE55,d18O_JARE55,'.')
% plot(sst_ncep_ACTIV,d18O_ACTIV,'.')
plot(sst_ncep_STRASSE,d18O_STRASSE,'.')
plot(sst_ncep_RARAAVIS,d18O_RARAAVIS,'.')
plot(sst_ncep_Bermuda,d18O_Bermuda,'.')

fig('units','inches','width',8,'height',6,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(sst_era_JARE55,d18O_JARE55,'.')
plot(sst_era_ACTIV,d18O_ACTIV,'.')
plot(sst_era_STRASSE,d18O_STRASSE,'.')
plot(sst_era_RARAAVIS,d18O_RARAAVIS,'.')
plot(sst_era_Bermuda,d18O_Bermuda,'.')

%% comapre to other measurements


load('./data/Uemura_2008_vapor.txt');
% load('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/Uemura_2008_vapor.txt');
vapor08_T=Uemura_2008_vapor(:,4);
vapor08_SST=Uemura_2008_vapor(:,5);
vapor08_h=Uemura_2008_vapor(:,6);
vapor08_d18O=Uemura_2008_vapor(:,7);
vapor08_dD=Uemura_2008_vapor(:,8);
vapor08_dxs=Uemura_2008_vapor(:,9);
%calculate RHn for Uemura 08 data
vapor08_TK0=vapor08_T-273.15;
vapor08_SSTK=vapor08_SST-273.15;

    vapor08_e_s_skin0 = (1000^-1).*exp(54.842763 - 6763.22 ./ vapor08_TK0 - 4.21 .* log(vapor08_TK0) + 0.000367 .* vapor08_TK0 +...
    tanh(0.0415 .* (vapor08_TK0 - 218.8)) .*  (53.878 - 1331.22 ./vapor08_TK0 - 9.44523 .* log(vapor08_TK0) + 0.014025 .* vapor08_TK0)) ;%with T in [K] and ew in [kPa]
    vapor08_e_s_SST0 = (1000^-1).*exp(54.842763 - 6763.22 ./ vapor08_SSTK - 4.21 .* log(vapor08_SSTK) + 0.000367 .* vapor08_SSTK +...
    tanh(0.0415 .* (vapor08_SSTK - 218.8)) .*  (53.878 - 1331.22 ./ vapor08_SSTK - 9.44523 .* log(vapor08_SSTK) + 0.014025 .* vapor08_SSTK)) ;%with T in [K] and ew in [kPa
vapor08_hn=real((vapor08_h.*vapor08_e_s_skin0)./vapor08_e_s_SST0);%normalized realitive humidity, using T0 as skin temp

load('./data/Uemura_2010_vapor.txt');
% load('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/Uemura_2010_vapor.txt');
vapor10_T=Uemura_2010_vapor(:,3);
vapor10_SST=Uemura_2010_vapor(:,4);
vapor10_h=Uemura_2010_vapor(:,5);
vapor10_hn=Uemura_2010_vapor(:,6);
vapor10_d18O=Uemura_2010_vapor(:,8);%these are the same samples as '08. but they got enriched
vapor10_d17O=Uemura_2010_vapor(:,9);
vapor10_d17Oxs=Uemura_2010_vapor(:,10);

% figure;
% hold on
% plot(vapor10_SST,vapor10_h,'.')
% plot(vapor10_SST,vapor10_hn,'o')
% plot(vapor08_SST,vapor08_h,'x')
% plot(vapor08_SST,vapor08_hn,'*')

water_vapor_isotope_data_LJF=load('./data/water_vapor_isotope_data_LJF_compiled.txt');
% water_vapor_isotope_data_LJF=load('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/water_vapor_isotope_data_LJF_compiled.txt');
vapor_LJF_d18O=water_vapor_isotope_data_LJF(:,3);
vapor_LJF_dD=water_vapor_isotope_data_LJF(:,4);
vapor_LJF_RH=water_vapor_isotope_data_LJF(:,5);
vapor_LJF_SST=water_vapor_isotope_data_LJF(:,6);
vapor_LJF_T=water_vapor_isotope_data_LJF(:,7);

ind=find(~isnan(vapor_LJF_SST) & ~isnan(vapor_LJF_RH) & ~isnan(vapor_LJF_dD) & ~isnan(vapor_LJF_d18O));

%%

fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(vapor08_SST,vapor08_d18O,'k.')
% plot(vapor_LJF_SST,vapor_LJF_d18O,'.','Color',[0.8 0.8 0.8])
plot(vapor_LJF_SST(ind),vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])
plot(sst_ncep_JARE55,d18O_JARE55,'.')
% plot(sst_ncep_ACTIV,d18O_ACTIV,'.')
plot(sst_ncep_STRASSE,d18O_STRASSE,'.')
plot(sst_ncep_RARAAVIS,d18O_RARAAVIS,'.')
% plot(sst_ncep_Bermuda,d18O_Bermuda,'.')


fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(vapor08_SST,vapor08_dD,'k.')
% plot(vapor_LJF_SST,vapor_LJF_dD,'.','Color',[0.8 0.8 0.8])
plot(vapor_LJF_SST(ind),vapor_LJF_dD(ind),'.','Color',[0.8 0.8 0.8])
plot(sst_ncep_JARE55,dD_JARE55,'.')
% plot(sst_ncep_ACTIV,dD_ACTIV,'.')
plot(sst_ncep_STRASSE,dD_STRASSE,'.')
plot(sst_ncep_RARAAVIS,dD_RARAAVIS,'.')
% plot(sst_ncep_Bermuda,dD_Bermuda,'.')


fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(vapor08_SST,vapor08_dD-8.*vapor08_d18O,'k.')
plot(vapor_LJF_SST(ind),vapor_LJF_dD(ind)-8.*vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])
plot(sst_ncep_JARE55,dD_JARE55-8.*d18O_JARE55,'.')
% plot(sst_ncep_ACTIV,dD_ACTIV-8.*d18O_ACTIV,'.')
plot(sst_ncep_STRASSE,dD_STRASSE-8.*d18O_STRASSE,'.')
plot(sst_ncep_RARAAVIS,dD_RARAAVIS-8.*d18O_RARAAVIS,'.')
% plot(sst_ncep_Bermuda,dD_Bermuda-8.*d18O_Bermuda,'.')




fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(sst_ncep_JARE55,dD_JARE55-8.*d18O_JARE55,'.')
plot(sst_era_JARE55,dD_JARE55-8.*d18O_JARE55,'.')

plot(vapor08_SST,vapor08_dD-8.*vapor08_d18O,'k.')
plot(vapor_LJF_SST(ind),vapor_LJF_dD(ind)-8.*vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])


fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(vapor08_h,vapor08_dD-8.*vapor08_d18O,'k.')
plot(vapor_LJF_RH(ind),vapor_LJF_dD(ind)-8.*vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])
plot(RH_JARE55.*100,dD_JARE55-8.*d18O_JARE55,'.')
%%
%~~~~~~~~~~~~~~~~~`

fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
subplot(131)
hold on
plot(Ta_JARE55,d18O_JARE55,'.')
plot(vapor_LJF_T(ind),vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])

plot(Ta_STRASSE,d18O_STRASSE,'.')
plot(Ta_RARAAVIS,d18O_RARAAVIS,'.')

scatter(vapor08_T,vapor08_d18O,'k','filled')


subplot(132)
hold on
plot(Ta_JARE55,dD_JARE55,'.')
plot(vapor_LJF_T(ind),vapor_LJF_dD(ind),'.','Color',[0.8 0.8 0.8])

plot(Ta_STRASSE,dD_STRASSE,'.')
plot(Ta_RARAAVIS,dD_RARAAVIS,'.')


scatter(vapor08_T,vapor08_dD,'k','filled')




subplot(133)
hold on

plot(Ta_JARE55,dD_JARE55-8.*d18O_JARE55,'.')
plot(vapor_LJF_T(ind),vapor_LJF_dD(ind)-8.*vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])


plot(Ta_STRASSE,dD_STRASSE-8.*d18O_STRASSE,'.')
plot(Ta_RARAAVIS,dD_RARAAVIS-8.*d18O_RARAAVIS,'.')


scatter(vapor08_T,vapor08_dD-8.*vapor08_d18O,'k','filled')


%~~~~~~~~~~~~~~~~~`
%%
fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on
plot(Ta_JARE55,dD_JARE55,'.')

fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on

plot(Ta_JARE55,d18O_JARE55,'.')
plot(vapor_LJF_T(ind),vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])
scatter(vapor08_T,vapor08_d18O,'k','filled')


fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on

plot(sst_ncep_JARE55,d18O_JARE55,'.')
plot(sst_era_alt_JARE55,d18O_JARE55,'.')

plot(vapor_LJF_SST(ind),vapor_LJF_d18O(ind),'.','Color',[0.8 0.8 0.8])
scatter(vapor08_SST,vapor08_d18O,'k','filled')


% 
% Ta_interp=-5:.1:25;
% [~, i_sort] = sort(Ta_JARE55);
% dD_JARE55_interp=interp1(Ta_JARE55(i_sort),dD_JARE55(i_sort),Ta_interp);
% d18O_JARE55_interp=interp1(Ta_JARE55(i_sort),d18O_JARE55(i_sort),Ta_interp);
% 
% fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
% hold on
% plot(Ta_JARE55,dD_JARE55,'.')
% plot(Ta_interp,dD_JARE55_interp)
% 
% fig('units','inches','width',11,'height',5,'font','Helvetica','fontsize',16,'border','on');
% hold on
% plot(Ta_JARE55,d18O_JARE55,'.')
% plot(Ta_interp,d18O_JARE55_interp)

Ta_binned=-5:.5:25;
window=.25;

for i =1:length(Ta_binned)
    dD_JARE55_binned_mean(i)=nanmean(dD_JARE55((Ta_JARE55>=(Ta_binned(i) -window)) & (Ta_JARE55<=(Ta_binned(i) +window)))); 
    dD_JARE55_binned_median(i)=nanmedian(dD_JARE55((Ta_JARE55>=(Ta_binned(i) -window)) & (Ta_JARE55<=(Ta_binned(i) +window)))); 
    d18O_JARE55_binned_mean(i)=nanmean(d18O_JARE55((Ta_JARE55>=(Ta_binned(i) -window)) & (Ta_JARE55<=(Ta_binned(i) +window)))); 
    d18O_JARE55_binned_median(i)=nanmedian(d18O_JARE55((Ta_JARE55>=(Ta_binned(i) -window)) & (Ta_JARE55<=(Ta_binned(i) +window)))); 
end
 
figure;
hold on
plot(Ta_JARE55,dD_JARE55,'.')
plot(Ta_binned,dD_JARE55_binned_mean)
plot(Ta_binned,dD_JARE55_binned_median)


 
figure;
hold on
plot(Ta_JARE55,d18O_JARE55,'.')
plot(Ta_binned,d18O_JARE55_binned_mean)
plot(Ta_binned,d18O_JARE55_binned_median)

% make means for each year
YEAR=2013:0.01:2016;
for i=1:length(YEAR)
idx=find(round(Ship_JARE55(2:end,1)./1e8,2)==YEAR(i) & diff(Ship_JARE55(:,1))<1e4);
d18O_JARE55_yearly(i)=nanmean(d18O_JARE55(idx));
dD_JARE55_yearly(i)=nanmean(dD_JARE55(idx));
Ta_JARE55_yearly(i)=nanmean(Ta_JARE55(idx));
RH_JARE55_yearly(i)=nanmean(RH_JARE55(idx));
clear idx
end


