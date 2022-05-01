function [fraction_i, fraction_l] = fraction_il_brm_H10(T,method)
% This is a function that calculates the fraction of ice and liquid water
% in clouds using polynomials from Hu et al 2010
% ref: Hu et al 2010, Occurrence, liquid water content, and fraction of supercooledwater clouds from combined CALIOP/IIR/MODIS measurements
% This can be compared to fraction_il_brm which uses an error function type
% fit following Kavanaugh and Cuffey
% T is the input temperature vector which should be in celcius
% Method can either be "mid", "iir", or "adj"
%% Housekeeping
if ~exist('method','var')
    method = 'mid';%default is the fit for "mid" level cloud temperature.
end
%%
%from ref: Hu et al 2010, Occurrence, liquid water content, and fraction of supercooledwater clouds from combined CALIOP/IIR/MODIS measurements

%midlayer temp
p_H10mid = 5.3608 + (0.4025.*T) + (0.08387.*T.^2) + (0.007182.*T.^3) + (2.39.*(10^-4).*T.^4) + (2.87.*(10^-6).*T.^5);
frac_l_H10mid=1./(1+exp(-p_H10mid));

%IIR temp
p_H10 = 5.2918 + (0.3694.*T) + (0.06635.*T.^2) + (0.006367.*T.^3) + (2.33.*(10^-4).*T.^4) + (2.97.*(10^-6).*T.^5);
frac_l_H10=1./(1+exp(-p_H10));

%adjusted, to look more like the "over land ice" IIR curve from Figure 6c
%Hu et al 2010
p_H10adj = 5.37 + (0.4025.*T) + (0.0847.*T.^2) + (0.007182.*T.^3) + (2.39.*(10^-4).*T.^4) + (2.87.*(10^-6).*T.^5);
frac_l_H10adj=1./(1+exp(-p_H10adj));

if method == 'mid'
    fraction_l = frac_l_H10mid;
    fraction_i = 1-fraction_l;
elseif method == 'iir'
    fraction_l = frac_l_H10;
    fraction_i = 1-fraction_l;
elseif method == 'adj'
    fraction_l = frac_l_H10adj;
    fraction_i = 1-fraction_l;
end



%  [fraction_i_KC03, fraction_l_KC03] = fraction_il_brm(T,0,-40);%used in kav and cuffey 2003
% fig('units','inches','width',6,'height',6,'font','Helvetica','fontsize',16,'border','on');
% hold on;
% plot(T,frac_l_H10,'k','Linewidth',2)
% plot(T,frac_l_H10mid,'--k','Linewidth',2)
% plot(T,frac_l_test,'Linewidth',2)
% plot(T,fraction_l_KC03,'Linewidth',2)