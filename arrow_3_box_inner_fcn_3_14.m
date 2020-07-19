    %'Arrow-3-box' function developed by Lorraine Lisiecki and Rachel Spratt, Geology
    %department, UCSB,2018, in order to simulate atmoshpheric fractionation of
    %the heavy isotope 18O as it is transported from the low to the high
    %latitudes via evaporation and subsequent Rayleigh fractionation
 
function  [ dX1,dX2,dX3]=arrow_3_box_inner_fcn_3_14(X1,X2,X3,n,k,overturning,mix_lo,mix_hi,evap_flux, horizontal,fraction,k2,lo_hi_expt,runoff);
    %input concentrations for box 1, box 2, box 3, iterations of the loop,
    %alpha1, overturning, low-latitude mixing, high-latitude
    %mixing,evaporative flux, horizontal mixing, fraction of transport
    %across 50 North, alpha2
    %
    %output change in concentration at each timestep for box 1, box 2, 
    %box 3, and atmospheric transport (same as box 1)
    
    %runoff=0.18;%Sverdrup strength of runoff
    %lo_hi_expt=0%.115;%Transport thru Bering Strait

    overturning_1_2=overturning-(fraction*evap_flux)-runoff-lo_hi_expt;%overturning flux values from lo latitude box to high latitude box in the ocean. 
    %We subtract the flux from the low to the high latitude box through the
    %atmosphere to  prevent the creation of mass
    standard=.1107/55.55;
    Y=((-20/1000)+1)*standard*55.55;  %runoff into high latitude surface box is -20 per mil
    Z=((-2/1000)+1)*standard*55.55; %transport thru Bering Strait is -2 per mil
    
        SPYR=86400*365;               % seconds per year

   atmospheric_transport=fraction*evap_flux*(k)*(fraction^(k2-1));
    V=[1.4858e+16 1.8745e+15  2.7476e+17]; %volumes of the 3 ocean boxes in m^3: in order; low, high, deep

    % for i=1:1:n; %number of time iterations for 18O concentration development through each cean box

    % Flux matrix
    F=[0 (overturning_1_2+horizontal+atmospheric_transport) mix_lo ; %row one, box one flux values  1,1 (from box 1 to box 1) 1,2 (from box 1 to box 2 o_t)  1,3 (from box 1 to box 3 o_th)
        horizontal       0             (mix_hi+overturning);        %row two, box two flux values  2,1 (from box 2 to box 1 t_o) 2,2 (from box 2 to box 2)  2,3 (from box 2 to box 3 t_th)
        (mix_lo+overturning)  mix_hi                  0];      %row 3,   box 3   flux values  3,1 (from box 3 to box 1 th_o) 3,2 (from box 3 to box 2 th_t)  3,3 (from box 3 to box 1)

    F = F.*1e6;%Sverdrup conversion [T X]=ode15s('box',Trange,Xinit); [T X]=ode15s('box',Trange,Xinit);
    %(mmol / [m^3])*(# Sv* 1e6 m^3/ yr) * (Sv)*(# seconds per yr)*(1/ ocean vol [#m^3])
runoff=runoff*1e6;
lo_hi_expt=lo_hi_expt*1e6;
    dX1=(((F(3,1)*X3+F(2,1)*X2)-(F(1,2)+F(1,3))*X1)-runoff*Y-lo_hi_expt*Z)*SPYR/V(1); %box 1, surface low-lat
    dX2=(((F(1,2)*X1+F(3,2)*X3+runoff*(Y)+lo_hi_expt*Z))-(F(2,3)+F(2,1))*X2)*SPYR/V(2); %box 2, surface high latitude
    dX3=(((F(1,3)*X1+F(2,3)*X2)-(F(3,1)+F(3,2))*X3))*SPYR/V(3);% box 3, deep ocean
  



