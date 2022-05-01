function [T_out] = Ts_to_Tc_2020(T_in,flip)
% This function converts surface temperature (T_s) to condensation temperature (T_c),
% over antarctica.
% I use a simple linear fit of T_c and T_S from Masson-Delmotte et al 2008.
% This is based on ERA data (c.f. Fig 8), and is in reasonable agreement
% MAR (a mesoscale model) data. The slope of the linear fit (0.65) is also in
% agreement with the slope of T_s to the inversion temperature from a
% couple older studies (0.67), though there appears to be an offset (the
% inversion temp is slightly warmer than the condensation temp). MD et al
% do not report the equation of the fit in the paper, so I estimate the
% y-intercept from the graph as T_c =-5. See below for modification
% The equation is T_c =  0.65 * T_s -5.

%The function can be used in reverse. If flip == 0, it will calculate the
%T_c from the T_s. If flip == 1 however, it will calculate a T_s from an
%input T_c. Default is flip = 0;

% UPDATE: I made a modified fit, with a y-intercept of -10. This drops the
% line in fig 8.a of MD08 (the ERA fit) a bit, putting in more in the
% middle of the cloud of the MAR data (though does not accoundt for
% apparently non-linearity in that data, naturally). This also makes the
% fit to my isotope model better...
% I think this is reasonable. The upper limit of the data in both MD08 Fig
% 8.a and 8.b seems to come from the inversion temperature. Condensation
% temp in the models is consistently lower than that. In panel a its a bit
% lower, in panel b, which is from a higher resolution model, its lower
% still. Making the simple assumption that its linear, we can use the slope
% from panel a and lower the intercept to fit it through the middle of the
% cloud in b. 


%update 2020: BRM
% T_surf=(T_cond- -8.1609)./0.6889;
% This is the best fit of observational data. 


%% Housekeeping

if ~exist('flip','var')
    flip = 0;%default is to calculate T_c from T_s
end
%%
% a=0.75;
% b=7;
% a=0.7;
% b=6;

% a=0.64;
% b=9;
% b=11;
% b=13;

% a=0.61;
% b=13;

% a=0.67;
% b=11;
% b=10;
% b=13;


a=0.69;
b=8.2;

% a=0.65;
% b=13;

%test based on MERRA2 data
%  a=0.59;  b=13;
if flip == 0%below T_in = surface Temp, and T_out = condensation temp. 
    
% T_out =  0.65 .* T_in -5;%fit based on ERA
% T_out =  0.67 .* T_in + -2;%inversion temp Jouzel and Merlivat 84 (intercept estmated from figure 8, MD 08)
% T_out =  0.65 .* T_in-9;%modified fit based on MAR
% T_out =  0.67 .* T_in-10;%modified fit based on MAR
% 
% T_out =  0.73 .* T_in-7;%modified fit based on MAR, also tuned a bit.
% T_out =  0.75 .* T_in-5;%modified fit based on MAR, also tuned a bit. %%%%%% THIS IS THE ONE I HAVE BEEN UING may 25 2017
% 
% T_out =   0.825 .*T_in -2.5;%fit that split difference between ERA and 1:1, more in line with spread from MAR;
T_out =   a .*T_in -b;%fit that split difference between ERA and 1:1, more in line with spread from MAR;

%TWO methods for adjust condensation temps at warmer sites:

    %1)
% T_out(T_out<T_in)=T_in(T_out<T_in);%correct the condensation temp for warm sites
    %2) 
% for i = 1:length(T_in)
%     if T_out(i)<T_in(i)
% T_out(i)=mean([T_in(i) T_out(i)]);%correct the condensation temp for warm sites
%     end
% end
% clear i


else %below T_out = surface Temp, and T_in = condensation temp. 
    
    
%     T_out = (T_in + 5)./0.65;%fit based on ERA
%     T_out = (T_in + 9) ./0.65;%modif ied fit base on MAR
% 
%         T_out = (T_in +7)./  0.73;%modified fit based on MAR, also tuned a bit.
%         T_out = (T_in +5)./  0.75;%modified fit based on MAR, also tuned a bit. %%%%%% THIS IS THE ONE I HAVE BEEN UING may 25 2017
%     T_out = (T_in +2.5)./  0.825;%fit that split difference between ERA and 1:1;
    T_out = (T_in +b)./  a;%fit that split difference between ERA and 1:1;

    
%TWO methods for adjust condensation temps at warmer sites:    
   %1) 
%  T_out(T_in<T_out)=T_in(T_in<T_out);%correct the condensation temp for warm sites
    
    %2)
% for i = 1:length(T_in)
%     if T_in(i)<T_out(i)
% T_out(i)=mean([T_in(i) T_out(i)]);%correct the condensation temp for warm sites
%     end
% end
% clear i

end
