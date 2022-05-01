function [rh0, delta_rh0, sst0, delta_sst0, rhn0, deltarhn0] = T_RH_RHn_2020(T0, SH, reanalysis, fit_method, season, remake_fits)
% BRM May 2020; updated from earlier scripts
% This is a function that returns a sea surface temp (SST), relative
% humidity (RH) and normalized relative humidity (RHn), given a user
% specified surface air temperautre (T0)

% The user can specify either the southern hemisphere (SH=1), or northern
% hemispehre (SH = 0). The T to SST fits are different between hemispheres.
% The T to RH fits will use data form both hemispheres.

% The user can specify which reanalysis data to use, either NCEP or
% (monthly) ERA Interim

% some notes (05/19/2020: I don't have the Seasonal ERA included yet. I
% have some questions about how to calcualte RHn correct for ERA data, see:
% T_RH_RHN_ERA_alt.m. I haven't included all the premade polynomials yet. I
% don't really like the polynomials-they're not good fit for high T0 (or
% too low T0).

% Notes: 06/2020 my favorite fits now are the spline fits. They have the
% smoothest derivatives (the binned means have very noisey derivatives),
% and lest edge problems. 

%% Housekeeping
if ~exist('SH','var')
    SH = 1;%1 (default) is Southern Hemisphere, 0 is Northern Hemisphere. 2 is average of both.
end

if ~exist('fit_method','var')
    fit_method = 'spline';%'spline' (default) is cubic spline (with noise tolerance); 'poly' = polynomial fits, 'bin' = binned means
end

if ~exist('season','var')
    season = 'annual';%default is annual average (only option for NCEP), other options (for ERA) are "DJF", "MAM", "JJA", "SON"
end

if ~exist('reanalysis','var')
    reanalysis = 'ncep';%default is ncep, other option is "era"
end

if ~exist('remake_fits','var')
remake_fits=0;%This will load previously calculated fits. Much faster. Turn this to 1 if you want to remake the different fits.
end

if strcmp(reanalysis,'ncep')
    folder='./data/ncep_data/';
elseif strcmp(reanalysis,'era')
    folder='./data/era_data/';
end


ord = 9; %order of polynomial if used.
window = 1; %size of window for binned means, if used.

%% Now do the conditional calcualtions


if remake_fits == 0 % Calculate values from previously made fits
    if  strcmp(fit_method,'bin')% Use fits from binned means
        %% First Load Data
        if SH == 1 %Southern Hemisphere
            if strcmp(reanalysis,'ncep') % for NCEP renalysis
                load(fullfile(folder,'ncep_SH_binned_models.mat'))        
            elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                load(fullfile(folder,'era_SH_binned_models.mat'))
            end
        elseif SH == 0 %Northern Hemisphere
            if strcmp(reanalysis,'ncep') % for NCEP renalysis
                load(fullfile(folder,'ncep_NH_binned_models.mat'))               
            elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                load(fullfile(folder,'era_NH_binned_models.mat'))
            end
        elseif SH == 2 %Average of both hemispheres 
             if strcmp(reanalysis,'ncep') % for NCEP renalysis
                load(fullfile(folder,'ncep_all_binned_models.mat'))               
            elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                load(fullfile(folder,'era_all_binned_models.mat'))
            end          
        end
        
        % Now make calcualtions
        % first SSTs
         sst0=interp1(T_model,sst_model_smooth,T0,'linear','extrap');
         delta_sst0=interp1(T_model,delta_sst_smooth,T0,'linear','extrap');
%         sst0=interp1(T_model,sst_model_smooth,T0);
%         delta_sst0=interp1(T_model,delta_sst_smooth,T0);
        % now RH
        rh0=interp1(T_model,rh_model_smooth,T0);
        delta_rh0=interp1(T_model,delta_rh_smooth,T0);  
        % fix values beyond interpolation domain
        rh0(T0>28)=interp1(T_model,rh_model_smooth,28);delta_rh0(T0>28)=interp1(T_model,delta_rh_smooth,28);
        rh0(T0<-10)=interp1(T_model,rh_model_smooth,-10);delta_rh0(T0>-10)=interp1(T_model,delta_rh_smooth,-10);
                 
        
    elseif strcmp(fit_method,'poly') % calculate fits from polynomials. Not my favorite for outside domain.
            %% First Load Data
                if SH == 1 %Southern Hemisphere
                    if strcmp(reanalysis,'ncep') % for NCEP renalysis
                        load(fullfile(folder,'ncep_SH_polynomials.mat'))        
                    elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                        load(fullfile(folder,'era_SH_polynomials.mat'))
                    end
                elseif SH == 0 %Northern Hemisphere
                    if strcmp(reanalysis,'ncep') % for NCEP renalysis
                        load(fullfile(folder,'ncep_NH_polynomials.mat'))               
                    elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                        load(fullfile(folder,'era_NH_polynomials.mat'))
                    end
                elseif SH == 2 %Average of both hemispheres 
                     if strcmp(reanalysis,'ncep') % for NCEP renalysis
                        load(fullfile(folder,'ncep_all_polynomials.mat'))               
                    elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                        load(fullfile(folder,'era_all_polynomials.mat'))
                    end          
                end
                
            [sst0, delta_sst0] = polyval(p_sst,T0,S_sst,mu_sst);
            [rh0, delta_rh0] = polyval(p_rh,T0,S_rh,mu_rh);

                     
    elseif strcmp(fit_method,'spline')% calculate fits from cubic spline
            %% First Load Data
                if SH == 1 %Southern Hemisphere
                    if strcmp(reanalysis,'ncep') % for NCEP renalysis
                        load(fullfile(folder,'ncep_spline_model_SH.mat'))        
                    elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                        load(fullfile(folder,'era_spline_model_SH.mat'))
                    end
                elseif SH == 0 %Northern Hemisphere
                    if strcmp(reanalysis,'ncep') % for NCEP renalysis
                        load(fullfile(folder,'ncep_spline_model_NH.mat'))               
                    elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                        load(fullfile(folder,'era_spline_model_NH.mat'))
                    end
                elseif SH == 2 %Average of both hemispheres 
                     if strcmp(reanalysis,'ncep') % for NCEP renalysis
                        load(fullfile(folder,'ncep_spline_model_all.mat'))               
                    elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
                        load(fullfile(folder,'era_spline_model_all.mat'))
                    end          
                end
 % Now make calcualtions
        % first SSTs
        sst0 =fnval(sp_sst,T0);
        delta_sst0 = interp1(T_model_spline,delta_sst_smooth,T0);
        % then rh
        rh0 =fnval(sp_rh,T0);
        delta_rh0 = interp1(T_model_spline,delta_rh_smooth,T0);  
        

        
    end
          %% now calculate rhn
%             if strcmp(reanalysis,'ncep') % for NCEP renalysis      
                TK0=T0+273.15;
                SSTK0=sst0+273.15;

                e_s_skin0 = (1000^-1).*exp(54.842763 - 6763.22 ./ TK0 - 4.21 .* log(TK0) + 0.000367 .* TK0 +...
                            tanh(0.0415 .* (TK0 - 218.8)) .*  (53.878 - 1331.22 ./ TK0 - 9.44523 .* log(TK0) + 0.014025 .* TK0)) ;%with T in [K] and ew in [kPa]
                e_s_SST0 = (1000^-1).*exp(54.842763 - 6763.22 ./ SSTK0 - 4.21 .* log(SSTK0) + 0.000367 .* SSTK0 +...
                            tanh(0.0415 .* (SSTK0 - 218.8)) .*  (53.878 - 1331.22 ./ SSTK0 - 9.44523 .* log(SSTK0) + 0.014025 .* SSTK0)) ;%with T in [K] and ew in [kPa

                rhn0=real((rh0.*e_s_skin0)./e_s_SST0);%normalized realitive humidity
                RHn2=real(((rh0+delta_rh0).*e_s_skin0)./e_s_SST0);%normalized realitive humidity
                deltarhn0=(RHn2-rhn0)/100;
                rhn0=rhn0/100;
%             elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
%                     %Now calculate RHN using ERA method:
%                     a1 = 611.21; %Pa;
%                     a3 = 17.502;
%                     a4 = 32.19; %K
%                     T_0 = 273.16; %K.
%                      e_sat_T0= a1*exp(a3 .* ((T0 - T_0) ./ (T0 - a4)));
%                      e_sat_SST0 = a1*exp(a3 .* ((sst0 - T_0) ./ (sst0 - a4)));
%                      rhn0=real((rh0.*e_sat_T0)./e_sat_SST0);%normalized realitive humidity
%                      RHn2=real(((rh0+delta_rh0).*e_sat_T0)./e_sat_SST0);%normalized realitive humidity
%                      deltarhn0=(RHn2-rhn0)/100;
%                      rhn0=rhn0/100;
%             end

elseif remake_fits == 1 % remake various fits. need to load data.

    
    
    %% Remake all fits
        T_model=-15:0.1:28; %T_model grid for binned means, if used.
%% Remake NCEP FITS
    if strcmp(reanalysis,'ncep') % for NCEP renalysis
            
            %load NCEP surface (skin) Temperature data
            lat_T=ncread(fullfile(folder,'yycompos.KqF1cVsn29.nc'),'lat');
            lon_T=ncread(fullfile(folder,'yycompos.KqF1cVsn29.nc'),'lon');
            T=ncread(fullfile(folder,'yycompos.KqF1cVsn29.nc'),'air');

            %load NCEP RH data
            rh=ncread(fullfile(folder,'yycompos.qhwjUBIKkV.nc'),'rhum');

            %load NCEP SST data (NOAA OI SST data set)
            SST=ncread(fullfile(folder,'yycompos.wQ7zPQ7Pct.nc'),'sst');
            lat_SST=ncread(fullfile(folder,'yycompos.wQ7zPQ7Pct.nc'),'lat');
            lon_SST=ncread(fullfile(folder,'yycompos.wQ7zPQ7Pct.nc'),'lon');
            SST(SST<-2)=nan;

            %set up lat lon grid
            [lats,lons] = meshgrid(lat_T,lon_T);
            [lats_SST,lons_SST] = meshgrid(lat_SST,lon_SST);

            %interp SST data onto surface temp grid
            SST_interp = interp2(lats_SST,lons_SST,SST,lats,lons);

            %make mask for ocean, will be NaN over continent, 1 elsewhere 
            mask=ones(size(SST_interp));
            mask(isnan(SST_interp))=NaN;
            mask(~isnan(SST_interp))=1;
            

            if SH == 1 %Southern Hemisphere
                T_line=reshape(T(lats<0),1,[]);
                SST_line=reshape(SST_interp(lats<0),1,[]);
                xsst=T_line(~isnan(SST_line));
                ysst=SST_line(~isnan(SST_line));
                %  xrh=T(~isnan(mask)& lats<0); % turn these on if you want to look at rh differences between the hemispheres. I don't think its needed.
                %  yrh=rh(~isnan(mask)& lats<0);
            elseif SH ==0 %Northern Hemisphere
                T_line=reshape(T(lats>0),1,[]);
                SST_line=reshape(SST_interp(lats>0),1,[]);
                xsst=T_line(~isnan(SST_line));
                ysst=SST_line(~isnan(SST_line));
                % xrh=T(~isnan(mask)& lats>0); % turn these on if you want to look at rh differences between the hemispheres. I don't think its needed.
                % yrh=rh(~isnan(mask)& lats>0);
            elseif SH ==2 %both hemispheres
                T_line=reshape(T,1,[]);
                SST_line=reshape(SST_interp,1,[]);
                xsst=T_line(~isnan(SST_line));
                ysst=SST_line(~isnan(SST_line));
            end
              xrh=T(~isnan(mask));% I group rh together for both hemispheres. there isn't a large difference between hems for rh, unlike sst
              yrh=rh(~isnan(mask));
              

              
%%
            if strcmp(fit_method,'poly')% polynomial method
            [p_sst,S_sst,mu_sst] = polyfit(xsst,ysst,ord);
            [sst0, delta_sst0] = polyval(p_sst,T0,S_sst,mu_sst);
            [p_rh,S_rh,mu_rh] = polyfit(xrh,yrh,ord);
            [rh0, delta_rh0] = polyval(p_rh,T0,S_rh,mu_rh);
            
            elseif strcmp(fit_method,'spline') %cubic spline method
                %SSTs
                [~,I_1] = sort(xsst); 
                x_1=xsst(I_1);
                y_1=ysst(I_1);
                [sp_sst] = csaps(x_1,y_1,0.001); %cubic spline interpoloation. smoothing parameter chosen for best fit/smooth derivative
                sst0 =fnval(sp_sst,T0);
                sst_diff=fnval(sp_sst,x_1);
                delta1=y_1-sst_diff;
                %rh
                [~,I_2] = sort(xrh); 
                x_2=xrh(I_2);
                y_2=yrh(I_2);
                [sp_rh] = csaps(x_2,y_2,0.0005); %cubic spline interpoloation. smoothing parameter chosen for best fit/smooth derivative
                rh0 =fnval(sp_rh,T0);
                rh_diff=fnval(sp_rh,x_2);
                delta2=y_2-rh_diff;                
               %calculate uncertainty in fit             
                window_spline=7;
                T_model_spline=-20:0.1:28;
                delta_sst=nan(size(T_model_spline));
                delta_rh=nan(size(T_model_spline));
                for i = 1:length(T_model_spline)
                    %SST
                    delta_sst(i)= nanstd(delta1((xsst)>=(T_model_spline(i)-window_spline) & (xsst)<=(T_model_spline(i)+window_spline)));
                    %RH
                    delta_rh(i)= nanstd(delta2((xrh)>=(T_model_spline(i)-window_spline) & (xrh)<=(T_model_spline(i)+window_spline)));
                end
              
                delta_sst_smooth=smooth(delta_sst,11)';
                delta_sst0 = interp1(T_model_spline,delta_sst_smooth,T0);
                delta_rh_smooth=smooth(delta_rh,11)';
                delta_rh0 = interp1(T_model_spline,delta_rh_smooth,T0);
               
            
            elseif  strcmp(fit_method,'bin') % binned means method
                sst_model=nan(size(T_model));delta_sst=nan(size(T_model));
                rh_model=nan(size(T_model));delta_rh=nan(size(T_model));
                for i = 1:length(T_model)
                    %SST
%                     sst_model(i)= nanmean(ysst((xsst)>=(T_model(i)-window) & (xsst)<=(T_model(i)+window)));
                  sst_model(i)= nanmedian(ysst((xsst)>=(T_model(i)-window) & (xsst)<=(T_model(i)+window)));
                    delta_sst(i)= nanstd(ysst((xsst)>=(T_model(i)-window) & (xsst)<=(T_model(i)+window)));
                    %RH
%                     rh_model(i)= nanmean(yrh((xrh)>=(T_model(i)-window) & (xrh)<=(T_model(i)+window)));
                    rh_model(i)= nanmedian(yrh((xrh)>=(T_model(i)-window) & (xrh)<=(T_model(i)+window)));
                    delta_rh(i)= nanstd(yrh((xrh)>=(T_model(i)-window) & (xrh)<=(T_model(i)+window)));
                end
                sst_model_smooth=smooth(smooth(sst_model,51),31);delta_sst_smooth=smooth(delta_sst,151);
                rh_model_smooth=smooth(rh_model,151);delta_rh_smooth=smooth(delta_rh,151);

                        % Now make calcualtions
                % first SSTs
                 sst0=interp1(T_model,sst_model_smooth,T0,'linear','extrap');
                 delta_sst0=interp1(T_model,delta_sst_smooth,T0,'linear','extrap');
%                 sst0=interp1(T_model,sst_model_smooth,T0);
%                 delta_sst0=interp1(T_model,delta_sst_smooth,T0);
                % now RH
                rh0=interp1(T_model,rh_model_smooth,T0);
                delta_rh0=interp1(T_model,delta_rh_smooth,T0);  
                % fix values beyond interpolation domain
                rh0(T0>28)=interp1(T_model,rh_model_smooth,28);delta_rh0(T0>28)=interp1(T_model,delta_rh_smooth,28);
                rh0(T0<-10)=interp1(T_model,rh_model_smooth,-10);delta_rh0(T0>-10)=interp1(T_model,delta_rh_smooth,-10);
            end
            %%
                %now calculate rhn, should work for all methods
                TK0=T0+273.15;
            	SSTK0=sst0+273.15;

                 e_s_skin0 = (1000^-1).*exp(54.842763 - 6763.22 ./ TK0 - 4.21 .* log(TK0) + 0.000367 .* TK0 +...
                   tanh(0.0415 .* (TK0 - 218.8)) .*  (53.878 - 1331.22 ./ TK0 - 9.44523 .* log(TK0) + 0.014025 .* TK0)) ;%with T in [K] and ew in [kPa]
                 e_s_SST0 = (1000^-1).*exp(54.842763 - 6763.22 ./ SSTK0 - 4.21 .* log(SSTK0) + 0.000367 .* SSTK0 +...
                   tanh(0.0415 .* (SSTK0 - 218.8)) .*  (53.878 - 1331.22 ./ SSTK0 - 9.44523 .* log(SSTK0) + 0.014025 .* SSTK0)) ;%with T in [K] and ew in [kPa

               rhn0=real((rh0.*e_s_skin0)./e_s_SST0);%normalized realitive humidity
               RHn2=real(((rh0+delta_rh0).*e_s_skin0)./e_s_SST0);%normalized realitive humidity
               deltarhn0=(RHn2-rhn0)/100;
               rhn0=rhn0/100;
        
    
            
    elseif strcmp(reanalysis,'era') 
        %% Make fits for ERA data from scratch
        %% Load ERA data
        lat=ncread(fullfile(folder,'ERA_int_monthly_data.nc'),'latitude');
        lon=ncread(fullfile(folder,'ERA_int_monthly_data.nc'),'longitude');
        [lats,lons] = meshgrid(lat,lon);

        T_dew=ncread(fullfile(folder,'ERA_int_monthly_data.nc'),'d2m');
        T=ncread(fullfile(folder,'ERA_int_monthly_data.nc'),'t2m');
        SST=ncread(fullfile(folder,'ERA_int_monthly_data.nc'),'sst');
        SP=ncread(fullfile(folder,'ERA_int_monthly_data.nc'),'sp');
        T_skin=ncread(fullfile(folder,'ERA_int_monthly_data.nc'),'skt');

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

        T_zonal=squeeze(mean(T,1));
        RH_zonal=squeeze(mean(RH,1));

        mask=ones(size(SST));
        mask(isnan(SST))=NaN;
        mask(~isnan(SST))=1;

        T_mask=T.*mask;
        T_dew_mask=T_dew.*mask;
        T_skin_mask=T_skin.*mask;
        RH_mask=RH.*mask;

        SST_climo=squeeze(nanmean(SST,3));
        T_mask_climo=squeeze(nanmean(T_mask,3));
        T_dew_mask_climo=squeeze(nanmean(T_dew_mask,3));
        T_skin_mask_climo=squeeze(nanmean(T_skin_mask,3));
        RH_mask_climo=squeeze(nanmean(RH_mask,3));
        
        if SH == 1 % Southern Hemisphere
            T_mask_climo_line=reshape(T_mask_climo(lats<0),1,[]);%climatology
            SST_climo_line=reshape(SST_climo(lats<0),1,[]);%climatology
            %xrh=reshape(T_mask_climo(~isnan(mask(:,:,1)) & lats<0),1,[]);% uncomment to segreate rh by hemisphere
            %yrh=reshape(RH_mask_climo(~isnan(mask(:,:,1)) & lats<0),1,[]);%uncomment to segreate rh by hemisphere
        elseif SH == 0 %Northern Hemisphere
            T_mask_climo_line=reshape(T_mask_climo(lats>0),1,[]);%climatology
            SST_climo_line=reshape(SST_climo(lats>0),1,[]);%climatology
            %xrh=reshape(T_mask_climo(~isnan(mask(:,:,1)) & lats>0 ),1,[]);%uncomment to segreate rh by hemisphere
            %yrh=reshape(RH_mask_climo(~isnan(mask(:,:,1)) & lats>0),1,[]);%uncomment to segreate rh by hemisphere
        elseif SH == 2 %Both Hemispheres
            T_mask_climo_line=reshape(T_mask_climo,1,[]);%climatology
            SST_climo_line=reshape(SST_climo,1,[]);%climatology
        end
        
        xsst=T_mask_climo_line(~isnan(SST_climo_line))-273.15;%climatology sst; NEED TO CONVERT TO C FROM K
        ysst=SST_climo_line(~isnan(SST_climo_line))-273.15;%climatology sst; NEED TO CONVERT TO C FROM K
        
        % now RH. if you want to 
        xrh=reshape(T_mask_climo(~isnan(mask(:,:,1))),1,[])-273.15;%climatology rh
        yrh=reshape(RH_mask_climo(~isnan(mask(:,:,1))),1,[]);%climatology rh
        
        if strcmp(fit_method,'poly') %Polynomial fit method
            [p_sst,S_sst,mu_sst] = polyfit(xsst,ysst,ord);
            [sst0, delta_sst0] = polyval(p_sst,T0,S_sst,mu_sst);
            [p_rh,S_rh,mu_rh] = polyfit(xrh,yrh,ord);
            [rh0, delta_rh0] = polyval(p_rh,T0,S_rh,mu_rh);
            
        elseif strcmp(fit_method,'spline')%Cubic spline method 
            %SST
            %OK this is dumb, but because the interpolation I use below doesn't like non unique values (I think) I need to add random,
            %very small noise to the data.
            xsst=xsst+(0.000001.*randn(size(xsst)));
            [~,I_1] = sort(xsst);
            xsst_interp=-20:0.001:28;
            ysst_interp=interp1(xsst(I_1),ysst(I_1),xsst_interp,'nearest');
            x_1=xsst_interp;
            y_1=ysst_interp;
            [sp_sst] = csaps(x_1,y_1,0.001);%fit a cubic spline, smoothing parameter is tuned
            sst0 =fnval(sp_sst,T0);%calculate sst0
            sst_diff=fnval(sp_sst,x_1);
            delta1=y_1-sst_diff;
            %rh
            xrh=xrh+(0.000001.*randn(size(xrh)));
            [~,I_2] = sort(xrh); 
            xrh_interp=-20:0.001:28;
            yrh_interp=interp1(xrh(I_2),yrh(I_2),xrh_interp,'nearest');
            x_2=xrh_interp;
            y_2=yrh_interp;
            [sp_rh] = csaps(x_2,y_2,0.001);%fit a cubic spline, smoothing parameter is tuned
            rh0 =fnval(sp_rh,T0);% calculate rh0
            rh_diff=fnval(sp_rh,x_2);
            delta2=y_2-rh_diff;
            %calculate uncertainty in fit
                window_spline=7;
                T_model_spline=-20:0.1:28;
                delta_sst=nan(size(T_model_spline));
                delta_rh=nan(size(T_model_spline));
                for i = 1:length(T_model_spline)
                    %SST
                    delta_sst(i)= nanstd(delta1((xsst_interp)>=(T_model_spline(i)-window_spline) & (xsst_interp)<=(T_model_spline(i)+window_spline)));
                    %RH
                    delta_rh(i)= nanstd(delta2((xrh_interp)>=(T_model_spline(i)-window_spline) & (xrh_interp)<=(T_model_spline(i)+window_spline)));
                end    
                delta_sst_smooth=smooth(delta_sst,11)';
                delta_sst0 = interp1(T_model_spline,delta_sst_smooth,T0);
                delta_rh_smooth=smooth(delta_rh,11)';
                delta_rh0 = interp1(T_model_spline,delta_rh_smooth,T0);
            
        elseif strcmp(fit_method,'bin')  %binned mean method
                sst_model=nan(size(T_model));delta_sst=nan(size(T_model));
                for i = 1:length(T_model)
                    %SST
                    sst_model(i)= nanmean(ysst((xsst)>=(T_model(i)-window) & (xsst)<=(T_model(i)+window)));
%                   sst_model(i)= nanmedian(y((x)>=(T_model(i)-window) & (x)<=(T_model(i)+window)));
                    delta_sst(i)= nanstd(ysst((xsst)>=(T_model(i)-window) & (xsst)<=(T_model(i)+window)));
                    %RH
                    rh_model(i)= nanmean(yrh((xrh)>=(T_model(i)-window) & (xrh)<=(T_model(i)+window)));
                    delta_rh(i)= nanstd(yrh((xrh)>=(T_model(i)-window) & (xrh)<=(T_model(i)+window)));
                end
                sst_model_smooth=smooth(sst_model,31);delta_sst_smooth=smooth(delta_sst,31);
                rh_model_smooth=smooth(rh_model,31);delta_rh_smooth=smooth(delta_rh,31);
                            
                
                % first SSTs
                 sst0=interp1(T_model,sst_model_smooth,T0,'linear','extrap');
                 delta_sst0=interp1(T_model,delta_sst_smooth,T0,'linear','extrap');
                  % sst0=interp1(T_model,sst_model_smooth,T0);
                  % delta_sst0=interp1(T_model,delta_sst_smooth,T0);
                % now RH
                rh0=interp1(T_model,rh_model_smooth,T0);
                delta_rh0=interp1(T_model,delta_rh_smooth,T0);  
                % fix values beyond interpolation domain
                rh0(T0>28)=interp1(T_model,rh_model_smooth,28);delta_rh0(T0>28)=interp1(T_model,delta_rh_smooth,28);
                rh0(T0<-10)=interp1(T_model,rh_model_smooth,-10);delta_rh0(T0>-10)=interp1(T_model,delta_rh_smooth,-10);
        end
        
%         %Now calculate RHN using ERA method:
%          e_sat_T0= a1*exp(a3 .* ((T0 - T_0) ./ (T0 - a4)));
%          e_sat_SST0 = a1*exp(a3 .* ((sst0 - T_0) ./ (sst0 - a4)));
%          rhn0=real((rh0.*e_sat_T0)./e_sat_SST0);%normalized realitive humidity
%          RHn2=real(((rh0+delta_rh0).*e_sat_T0)./e_sat_SST0);%normalized realitive humidity
%          deltarhn0=(RHn2-rhn0)/100;
%          rhn0=rhn0/100;
% 
          %% now calculate rhn
%             if strcmp(reanalysis,'ncep') % for NCEP renalysis      
                TK0=T0+273.15;
                SSTK0=sst0+273.15;

                e_s_skin0 = (1000^-1).*exp(54.842763 - 6763.22 ./ TK0 - 4.21 .* log(TK0) + 0.000367 .* TK0 +...
                            tanh(0.0415 .* (TK0 - 218.8)) .*  (53.878 - 1331.22 ./ TK0 - 9.44523 .* log(TK0) + 0.014025 .* TK0)) ;%with T in [K] and ew in [kPa]
                e_s_SST0 = (1000^-1).*exp(54.842763 - 6763.22 ./ SSTK0 - 4.21 .* log(SSTK0) + 0.000367 .* SSTK0 +...
                            tanh(0.0415 .* (SSTK0 - 218.8)) .*  (53.878 - 1331.22 ./ SSTK0 - 9.44523 .* log(SSTK0) + 0.014025 .* SSTK0)) ;%with T in [K] and ew in [kPa

                rhn0=real((rh0.*e_s_skin0)./e_s_SST0);%normalized realitive humidity
                RHn2=real(((rh0+delta_rh0).*e_s_skin0)./e_s_SST0);%normalized realitive humidity
                deltarhn0=(RHn2-rhn0)/100;
                rhn0=rhn0/100;
%             elseif strcmp(reanalysis,'era') && strcmp(season,'annual') % for ERA renalysis
%                     %Now calculate RHN using ERA method:
%                     a1 = 611.21; %Pa;
%                     a3 = 17.502;
%                     a4 = 32.19; %K
%                     T_0 = 273.16; %K.
%                      e_sat_T0= a1*exp(a3 .* ((T0 - T_0) ./ (T0 - a4)));
%                      e_sat_SST0 = a1*exp(a3 .* ((sst0 - T_0) ./ (sst0 - a4)));
%                      rhn0=real((rh0.*e_sat_T0)./e_sat_SST0);%normalized realitive humidity
%                      RHn2=real(((rh0+delta_rh0).*e_sat_T0)./e_sat_SST0);%normalized realitive humidity
%                      deltarhn0=(RHn2-rhn0)/100;
%                      rhn0=rhn0/100;
%             end

end
end



    


        
        
