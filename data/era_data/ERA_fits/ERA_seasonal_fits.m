% ncdisp('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc')
%load ERA data
%data from http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/requests/netcdf/58bdfa3eff29d0668db5aa72/
lat=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc','latitude');
lon=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc','longitude');
[lats,lons] = meshgrid(lat,lon);

T_dew=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc','d2m');
T=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc','t2m');
SST=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc','sst');
SP=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc','sp');
T_skin=ncread('/Users/Bradley/Documents/work/wais_dxs_paper/WDC_dxs_2016/redefining_dxs_in_antarctica/papers_and_data/ERA_int_monthly_data.nc','skt');


%for calculating relative humidity see notes at http://www.ecmwf.int/en/does-era-40-dataset-contain-near-surface-humidity-data
%http://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf
a1 = 611.21; %Pa;
a3 = 17.502;
a4 = 32.19; %K
T_0 = 273.16; %K.
% e_sat = a1*exp(a3*((T ? T_0)/(T ? a4)));
e_sat_d = a1*exp(a3 .* ((T_dew - T_0) ./ (T_dew - a4)));
e_sat = a1*exp(a3 .* ((T - T_0) ./ (T - a4)));
% q_sat = ((R_dry/R_vap).*e_sat) ./ (SP - (1- (R_dry/R_vap).*e_sat));

RH=100 .* (e_sat_d ./ e_sat);

mask=ones(size(SST));
mask(isnan(SST))=NaN;
mask(~isnan(SST))=1;

T_mask=T.*mask;
T_dew_mask=T_dew.*mask;
T_skin_mask=T_skin.*mask;
RH_mask=RH.*mask;
% 
% T_mask_line=reshape(T_mask,1,[]);
% RH_mask_line=reshape(RH_mask,1,[]);
% SST_line=reshape(SST,1,[]);

% SST_climo=squeeze(nanmean(SST,3));
% T_mask_climo=squeeze(nanmean(T_mask,3));
% T_dew_mask_climo=squeeze(nanmean(T_dew_mask,3));
% T_skin_mask_climo=squeeze(nanmean(T_skin_mask,3));
% RH_mask_climo=squeeze(nanmean(RH_mask,3));


%%now make model for T to SST


T_model=-20:28; T_model=T_model+273.15;

%order of polynomial fits used throughout. ord = 9 seems to do well.
ord = 9;

%% look into seasonal changes?
T_DJF=[];
RH_DJF=[];
Tskin_DJF=[];
SST_DJF=[];

T_MAM=[];
RH_MAM=[];
Tskin_MAM=[];
SST_MAM=[];

T_JJA=[];
RH_JJA=[];
Tskin_JJA=[];
SST_JJA=[];


T_SON=[];
RH_SON=[];
Tskin_SON=[];
SST_SON=[];
%%

SH = 0;
%% SOUTHERN HEMISPHERE
if SH==1  
clear w x y z
%DJF
for i=1:12:361%January, for 30 years
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_DJF=[SST_DJF; w(~isnan(mask(:,:,1)) & lats<0)];
T_DJF=[T_DJF; x(~isnan(mask(:,:,1)) & lats<0)];
RH_DJF=[RH_DJF; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_DJF=[Tskin_DJF; z(~isnan(mask(:,:,1)) & lats<0)];
end
clear w x y z
for i=2:12:362%February
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_DJF=[SST_DJF; w(~isnan(mask(:,:,1)) & lats<0)];
T_DJF=[T_DJF; x(~isnan(mask(:,:,1)) & lats<0)];
RH_DJF=[RH_DJF; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_DJF=[Tskin_DJF; z(~isnan(mask(:,:,1)) & lats<0)];
end
clear w x y z
for i=12:12:372%December
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_DJF=[SST_DJF; w(~isnan(mask(:,:,1)) & lats<0)];
T_DJF=[T_DJF; x(~isnan(mask(:,:,1)) & lats<0)];
RH_DJF=[RH_DJF; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_DJF=[Tskin_DJF; z(~isnan(mask(:,:,1)) & lats<0)];
end


clear w x y z
%MAM
for i=3:12:363%March, for 30 years
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_MAM=[SST_MAM; w(~isnan(mask(:,:,1)) & lats<0)];
T_MAM=[T_MAM; x(~isnan(mask(:,:,1)) & lats<0)];
RH_MAM=[RH_MAM; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_MAM=[Tskin_MAM; z(~isnan(mask(:,:,1)) & lats<0)];
end
clear w x y z
for i=4:12:364%April
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_MAM=[SST_MAM; w(~isnan(mask(:,:,1)) & lats<0)];
T_MAM=[T_MAM; x(~isnan(mask(:,:,1)) & lats<0)];
RH_MAM=[RH_MAM; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_MAM=[Tskin_MAM; z(~isnan(mask(:,:,1)) & lats<0)];
end
clear w x y z
for i=5:12:365%May
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_MAM=[SST_MAM; w(~isnan(mask(:,:,1)) & lats<0)];
T_MAM=[T_MAM; x(~isnan(mask(:,:,1)) & lats<0)];
Tskin_MAM=[Tskin_MAM; z(~isnan(mask(:,:,1)) & lats<0)];
RH_MAM=[RH_MAM; y(~isnan(mask(:,:,1)) & lats<0)];
end

%JJA
clear w x y z
for i=6:12:366%June
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_JJA=[SST_JJA; w(~isnan(mask(:,:,1)) & lats<0)];
T_JJA=[T_JJA; x(~isnan(mask(:,:,1)) & lats<0)];
RH_JJA=[RH_JJA; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_JJA=[Tskin_JJA; z(~isnan(mask(:,:,1)) & lats<0)];
end

clear w x y z
for i=7:12:367%July
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_JJA=[SST_JJA; w(~isnan(mask(:,:,1)) & lats<0)];
T_JJA=[T_JJA; x(~isnan(mask(:,:,1)) & lats<0)];
RH_JJA=[RH_JJA; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_JJA=[Tskin_JJA; z(~isnan(mask(:,:,1)) & lats<0)];
end
clear w x y z
for i=8:12:368%August
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_JJA=[SST_JJA; w(~isnan(mask(:,:,1)) & lats<0)];
T_JJA=[T_JJA; x(~isnan(mask(:,:,1)) & lats<0)];
RH_JJA=[RH_JJA; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_JJA=[Tskin_JJA; z(~isnan(mask(:,:,1)) & lats<0)];
end


%SON
clear w x y z
for i=9:12:369%September
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_SON=[SST_SON; w(~isnan(mask(:,:,1)) & lats<0)];
T_SON=[T_SON; x(~isnan(mask(:,:,1)) & lats<0)];
RH_SON=[RH_SON; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_SON=[Tskin_SON; z(~isnan(mask(:,:,1)) & lats<0)];
end
clear w x y z
for i=10:12:370%October
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_SON=[SST_SON; w(~isnan(mask(:,:,1)) & lats<0)];
T_SON=[T_SON; x(~isnan(mask(:,:,1)) & lats<0)];
RH_SON=[RH_SON; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_SON=[Tskin_SON; z(~isnan(mask(:,:,1)) & lats<0)];
end
clear w x y z
for i=11:12:371%November
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_SON=[SST_SON; w(~isnan(mask(:,:,1)) & lats<0)];
T_SON=[T_SON; x(~isnan(mask(:,:,1)) & lats<0)];
RH_SON=[RH_SON; y(~isnan(mask(:,:,1)) & lats<0)];
Tskin_SON=[Tskin_SON; z(~isnan(mask(:,:,1)) & lats<0)];
end

% NOTHERN HEMISPHERE
elseif SH == 0
clear w x y z
%DJF
for i=1:12:361%January, for 30 years
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_DJF=[SST_DJF; w(~isnan(mask(:,:,1)) & lats>0)];
T_DJF=[T_DJF; x(~isnan(mask(:,:,1)) & lats>0)];
RH_DJF=[RH_DJF; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_DJF=[Tskin_DJF; z(~isnan(mask(:,:,1)) & lats>0)];
end
clear w x y z
for i=2:12:362%February
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_DJF=[SST_DJF; w(~isnan(mask(:,:,1)) & lats>0)];
T_DJF=[T_DJF; x(~isnan(mask(:,:,1)) & lats>0)];
RH_DJF=[RH_DJF; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_DJF=[Tskin_DJF; z(~isnan(mask(:,:,1)) & lats>0)];
end
clear w x y z
for i=12:12:372%December
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_DJF=[SST_DJF; w(~isnan(mask(:,:,1)) & lats>0)];
T_DJF=[T_DJF; x(~isnan(mask(:,:,1)) & lats>0)];
RH_DJF=[RH_DJF; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_DJF=[Tskin_DJF; z(~isnan(mask(:,:,1)) & lats>0)];
end


clear w x y z
%MAM
for i=3:12:363%March, for 30 years
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_MAM=[SST_MAM; w(~isnan(mask(:,:,1)) & lats>0)];
T_MAM=[T_MAM; x(~isnan(mask(:,:,1)) & lats>0)];
RH_MAM=[RH_MAM; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_MAM=[Tskin_MAM; z(~isnan(mask(:,:,1)) & lats>0)];
end
clear w x y z
for i=4:12:364%April
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_MAM=[SST_MAM; w(~isnan(mask(:,:,1)) & lats>0)];
T_MAM=[T_MAM; x(~isnan(mask(:,:,1)) & lats>0)];
RH_MAM=[RH_MAM; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_MAM=[Tskin_MAM; z(~isnan(mask(:,:,1)) & lats>0)];
end
clear w x y z
for i=5:12:365%May
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_MAM=[SST_MAM; w(~isnan(mask(:,:,1)) & lats>0)];
T_MAM=[T_MAM; x(~isnan(mask(:,:,1)) & lats>0)];
Tskin_MAM=[Tskin_MAM; z(~isnan(mask(:,:,1)) & lats>0)];
RH_MAM=[RH_MAM; y(~isnan(mask(:,:,1)) & lats>0)];
end

%JJA
clear w x y z
for i=6:12:366%June
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_JJA=[SST_JJA; w(~isnan(mask(:,:,1)) & lats>0)];
T_JJA=[T_JJA; x(~isnan(mask(:,:,1)) & lats>0)];
RH_JJA=[RH_JJA; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_JJA=[Tskin_JJA; z(~isnan(mask(:,:,1)) & lats>0)];
end

clear w x y z
for i=7:12:367%July
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_JJA=[SST_JJA; w(~isnan(mask(:,:,1)) & lats>0)];
T_JJA=[T_JJA; x(~isnan(mask(:,:,1)) & lats>0)];
RH_JJA=[RH_JJA; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_JJA=[Tskin_JJA; z(~isnan(mask(:,:,1)) & lats>0)];
end
clear w x y z
for i=8:12:368%August
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_JJA=[SST_JJA; w(~isnan(mask(:,:,1)) & lats>0)];
T_JJA=[T_JJA; x(~isnan(mask(:,:,1)) & lats>0)];
RH_JJA=[RH_JJA; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_JJA=[Tskin_JJA; z(~isnan(mask(:,:,1)) & lats>0)];
end


%SON
clear w x y z
for i=9:12:369%September
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_SON=[SST_SON; w(~isnan(mask(:,:,1)) & lats>0)];
T_SON=[T_SON; x(~isnan(mask(:,:,1)) & lats>0)];
RH_SON=[RH_SON; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_SON=[Tskin_SON; z(~isnan(mask(:,:,1)) & lats>0)];
end
clear w x y z
for i=10:12:370%October
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_SON=[SST_SON; w(~isnan(mask(:,:,1)) & lats>0)];
T_SON=[T_SON; x(~isnan(mask(:,:,1)) & lats>0)];
RH_SON=[RH_SON; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_SON=[Tskin_SON; z(~isnan(mask(:,:,1)) & lats>0)];
end
clear w x y z
for i=11:12:371%November
ind=i;
w=squeeze(SST(:,:,ind));
x=squeeze(T_mask(:,:,ind));
y=squeeze(RH_mask(:,:,ind));
z=squeeze(T_skin_mask(:,:,ind));
SST_SON=[SST_SON; w(~isnan(mask(:,:,1)) & lats>0)];
T_SON=[T_SON; x(~isnan(mask(:,:,1)) & lats>0)];
RH_SON=[RH_SON; y(~isnan(mask(:,:,1)) & lats>0)];
Tskin_SON=[Tskin_SON; z(~isnan(mask(:,:,1)) & lats>0)];
end
end

%%
[p_rh_DJF,S_rh_DJF,mu_rh_DJF] = polyfit(T_DJF,RH_DJF,ord);
[p_Tskin_DJF,S_Tskin_DJF,mu_Tskin_DJF] = polyfit(T_DJF,Tskin_DJF,ord);
[p_sst_DJF,S_sst_DJF,mu_sst_DJF] = polyfit(T_DJF,SST_DJF,ord);

[p_rh_MAM,S_rh_MAM,mu_rh_MAM] = polyfit(T_MAM,RH_MAM,ord);
[p_Tskin_MAM,S_Tskin_MAM,mu_Tskin_MAM] = polyfit(T_MAM,Tskin_MAM,ord);
[p_sst_MAM,S_sst_MAM,mu_sst_MAM] = polyfit(T_MAM,SST_MAM,ord);

[p_rh_JJA,S_rh_JJA,mu_rh_JJA] = polyfit(T_JJA,RH_JJA,ord);
[p_Tskin_JJA,S_Tskin_JJA,mu_Tskin_JJA] = polyfit(T_JJA,Tskin_JJA,ord);
[p_sst_JJA,S_sst_JJA,mu_sst_JJA] = polyfit(T_JJA,SST_JJA,ord);

[p_rh_SON,S_rh_SON,mu_rh_SON] = polyfit(T_SON,RH_SON,ord);
[p_Tskin_SON,S_Tskin_SON,mu_Tskin_SON] = polyfit(T_SON,Tskin_SON,ord);
[p_sst_SON,S_sst_SON,mu_sst_SON] = polyfit(T_SON,SST_SON,ord);


[rh_DJF_model, delta_rh_DJF] = polyval(p_rh_DJF,T_model(13:end),S_rh_DJF,mu_rh_DJF);
[rh_MAM_model, delta_rh_MAM] = polyval(p_rh_MAM,T_model,S_rh_MAM,mu_rh_MAM);
[rh_JJA_model, delta_rh_JJA] = polyval(p_rh_JJA,T_model,S_rh_JJA,mu_rh_JJA);
[rh_SON_model, delta_rh_SON] = polyval(p_rh_SON,T_model,S_rh_SON,mu_rh_SON);



figure;hold on;
plot(T_JJA,RH_JJA,'b.');
plot(T_DJF,RH_DJF,'r.');
plot(T_MAM,RH_MAM,'c.');
plot(T_SON,RH_SON,'k.');

plot(T_model(13:end),rh_DJF_model);
plot(T_model,rh_JJA_model);
plot(T_model,rh_MAM_model);
plot(T_model,rh_SON_model);
% plot(T_model,rh_model,'k');

if SH == 1
save('ERA_fits_SH_DJF.mat','p_sst_DJF','S_sst_DJF','mu_sst_DJF','p_Tskin_DJF','S_Tskin_DJF','mu_Tskin_DJF','p_rh_DJF','S_rh_DJF','mu_rh_DJF');
save('ERA_fits_SH_MAM.mat','p_sst_MAM','S_sst_MAM','mu_sst_MAM','p_Tskin_MAM','S_Tskin_MAM','mu_Tskin_MAM','p_rh_MAM','S_rh_MAM','mu_rh_MAM');
save('ERA_fits_SH_JJA.mat','p_sst_JJA','S_sst_JJA','mu_sst_JJA','p_Tskin_JJA','S_Tskin_JJA','mu_Tskin_JJA','p_rh_JJA','S_rh_JJA','mu_rh_JJA');
save('ERA_fits_SH_SON.mat','p_sst_SON','S_sst_SON','mu_sst_SON','p_Tskin_SON','S_Tskin_SON','mu_Tskin_SON','p_rh_SON','S_rh_SON','mu_rh_SON');
elseif SH ==0
save('ERA_fits_NH_DJF.mat','p_sst_DJF','S_sst_DJF','mu_sst_DJF','p_Tskin_DJF','S_Tskin_DJF','mu_Tskin_DJF','p_rh_DJF','S_rh_DJF','mu_rh_DJF');
save('ERA_fits_NH_MAM.mat','p_sst_MAM','S_sst_MAM','mu_sst_MAM','p_Tskin_MAM','S_Tskin_MAM','mu_Tskin_MAM','p_rh_MAM','S_rh_MAM','mu_rh_MAM');
save('ERA_fits_NH_JJA.mat','p_sst_JJA','S_sst_JJA','mu_sst_JJA','p_Tskin_JJA','S_Tskin_JJA','mu_Tskin_JJA','p_rh_JJA','S_rh_JJA','mu_rh_JJA');
save('ERA_fits_NH_SON.mat','p_sst_SON','S_sst_SON','mu_sst_SON','p_Tskin_SON','S_Tskin_SON','mu_Tskin_SON','p_rh_SON','S_rh_SON','mu_rh_SON');
end

