
clear all
    close all
    tic% Calculate the length of the run
    format short
    conversion1 = 0.110712368610173;%Mean 18O in the Atlantic

    %%     traces the heavy isotope for each timestep
    
  lo_hi_expt=.115;
  lo_hi_expt =[0.092 lo_hi_expt 0.138];
  
    
    runoff=0.18;
    runoff=[0.045 0.09 0.18];

    start_T=288;% Alpha 1 as a function of low latitude sea surface temperature
    end_T=290 ;

    T_1=[283 287 290]';

    alpha=(1./((exp((1137./(T_1.^2))-(.4156./(T_1))-.0020667))))' ;%(Majoube, 1971)
    d=length(alpha);% d gives us the same length for all experiments

    %start=5; %base overturning before loop
    %strongest=8;
    ov=[3 5 5.6];%start:.5:strongest;  %overturning flux value in Sverdrups: values are for each loop iteration 1-3
   ov=ov';
    e=length(ov);

    mix_l_start=2.5; %low-latitude exchnge with deep ocean
    mix_low_latitude=[1.75 3.5  7];
    f=length(mix_low_latitude);

    mix_h_start=5; %high latitude exchange with the deep ocean
    mix_high_latitude=[3.5  6.5  10];
    g=length(mix_high_latitude);

    evap=[0.72*.31   0.82*.31 .31];%20:.03:.30;%Evaporative flux over the low latitude surface ocean
    h=length(evap);
    %fractn_1=.36%round(.14/evap(end),2)
%fractn_end=.39%round(.14/evap(1),2)

    mix_horiz=[.45 .5 .55];% horizontal exchange between the two surface boxes (Eckman transport)
    l= length(mix_horiz);

   
    fract=[.19 .29 .39];%fractn_1:.03:fractn_end; % fraction of transport calculated between 40-50 degrees north and transported across 50 degrees North to the high latitudes
    p=length(fract);

    starting_T=start_T-3; %Alpha 2 as a function of temperature
    ending_T=end_T-3;
    T_2=(T_1-3)';

    alph2=(exp((1137./(T_2.^2))-(.4156./(T_2))-.0020667))';%(Majoube, 1971)

    v=length(alph2);
    ii=length(lo_hi_expt);
    ix=length(runoff);

    for q=1:e;%overturn from 13.5-17.5

        for r=1:f;%mix low latitude 4

            for s=1:g;%mix high latitude 5
                for t=1:d; %next change in alpha1 from ~0.99 to ~0.9908
                    for u=1:h; %evaporative flux over low latitude surface
                        for m=1:l;%mix horizontal
                            for o=1:p;%fraction sent to high latitudes ~40%
                                for w=1:v; %alpha 2 ~1.014
                                    for x=1:ii;
                                        for y=1:ix;
                                        
                                    n=900;% timestep in years (approximation)

                                    overturning=ov(q);
                                    k=alpha(t); %fractionation value for heavy isotope in low latitudes
                                    mix_lo=mix_low_latitude(r);
                                    mix_hi=mix_high_latitude(s);
                                    evap_flux=evap(u);
                                    horizontal=mix_horiz(m);
                                    fraction=fract(o);%~Rayleigh fractionation value for heavy isotope in high latitudes
                                    k2=alph2(w);
                                    l_h=lo_hi_expt(x);
                                    run=runoff(y);

                                    standard=.1107/55.55;   %Standard mean ocean water value in mmol/m^3

                                    dH1=zeros(n,1);
                                    dH2=zeros(n,1);
                                    dH3=zeros(n,1);
                                    dHf=zeros(n,1);
                                    H1=zeros(n,1);
                                    H2=zeros(n,1);
                                    H3=zeros(n,1);
                                    Hf=zeros(n,1);
                                    %

                                    H1(1)=  conversion1 ;  %Initial 18O concentration in mmol/m^3 for low latitude ocean box
                                    H2(1)=  conversion1 ;  %Initial 18O concentration in mmol/m^3 for high latitude ocean box
                                    H3(1)=  conversion1 ;  %Initial 18O concentration in mmol/m^3 for deep ocean box
                                    Hf(1)=H1(1)*k*fraction^(k2-1); %Calculating initial atmospheric flux value


                                    heavy_isotope=@arrow_3_box_inner_fcn_3_14;
                                    % Creates an alias called 'heavy isotope'-for
                                    % the box function
                                    d18O1=zeros(3,3,3,3,3,3,3,3,3,3);
                                    d18O2=zeros(3,3,3,3,3,3,3,3,3,3);
                                    d18O3=zeros(3,3,3,3,3,3,3,3,3,3);
                                    d18Of=zeros(3,3,3,3,3,3,3,3,3,3);
                                    for i=1:n;
                                        [dH1(i),dH2(i),dH3(i)]=heavy_isotope(H1(i),H2(i),H3(i),n,k,overturning,mix_lo,mix_hi,evap_flux, horizontal,fraction,k2,l_h,run);

                                        %Calls the box function with the
                                        %initial 18O values for each box
                                        %returns the '-change in'- 18O
                                        %for each timestep for each box



                                        %This little loop generates the
                                        %isotope value at the next timestep

                                        H1(i+1)=H1(i)+ dH1(i);%Box 1 evolving with each timestep
                                        H2(i+1)=H2(i)+ dH2(i);%Box 2 evolving with each timestep
                                        H3(i+1)=H3(i)+ dH3(i);%Box 3 evolving with each timestep
                                     Hf(i)=H1(i)*k*(fraction^(k2-1));%This generates the atmospheric transport value in per mil

                                   






                                    %%     traces the light isotope for each timestep


                                    box1=(((H1(i+1)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box1
                                    box2=(((H2(i+1)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box2
                                    box3=(((H3(i+1)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box3
                                    boxf=(((Hf(i)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box3
one(q,r,s,t,u,m,o,w,x,y)= box1(end);
two(q,r,s,t,u,m,o,w,x,y)= box2(end);
three(q,r,s,t,u,m,o,w,x,y)= box3(end);
fb(q,r,s,t,u,m,o,w,x,y)= boxf(end);



                                    end
                               
                                        end
                             
                        end
                    end
                end
            end

        end

            end
        end
        end
    
    
      % d18O1(q,r,s,t,u,m,o,w,x,y)= H1(end)
    d18O1(q,r,s,t,u,m,o,w,x,y)=(((H1(end)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box1
    d18O2(q,r,s,t,u,m,o,w,x,y)=(((H2(end)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box2
    d18O3(q,r,s,t,u,m,o,w,x,y)=(((H3(end)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box3
    d18Of(q,r,s,t,u,m,o,w,x,y)=(((Hf(end)/55.55)./standard)-1)*1000;%d18O values for overturn experiment box3



    %%


    figure (1)
    plot(((H1/55.55)/standard-1)*1000,'r','linewidth',3);
    hold on
    plot(((H2/55.55)/standard-1)*1000,'b','linewidth',3);
    hold on
    plot(((H3/55.55)/standard-1)*1000,'k','linewidth',3);
    %hold on
    %plot(((Hf/55.55)/standard-1)*1000,'o');
    title('single box evolution in years to steady state');
    ylabel('delta 18O');
    legend('low','high','deep');%,'atm fr')


    figure (2)
    plot(((H1/55.55)/standard-1)*1000,'r','linewidth',3);
    hold on
    plot(((H2/55.55)/standard-1)*1000,'b','linewidth',3);
    hold on
    plot(((H3/55.55)/standard-1)*1000,'k','linewidth',3);
    hold on
    plot(((Hf/55.55)/standard-1)*1000,'o');
    title('single box evolution in years to steady state');
    ylabel('delta 18O');
    xlabel('years to steady state');
    legend('low','high','deep');%,'atm fr')


    



    %%

    end
    figure (4)% low



  
    plot(box1(1,:,1,1,1,1),'k','linewidth',2);
    hold on
    plot( (box2(1,:,1,1,1,1)),'y','linewidth',2);
    hold on
    plot( (box3(1,:,1,1,1,1)),'g','linewidth',2);
    hold on

    xlabel('low latitude mixing');
    title('increasing low latitude mixing');

    legend('low','high','deep');

    % hold on
    toc
    
  ov=    [ov(1)                 ov(1)                ov(1)                       ov(1)                  ov(3)                 ov(1)          ov(1)                    ov(1)                ov(3)               ov(2)              ov(3)              ov(1)                     ov(3)                ov(3)                   ov(1)                       ov(1)                ov(3)]';
  l_l=  [mix_low_latitude(1)  mix_low_latitude(1)   mix_low_latitude(1)   mix_low_latitude(1) mix_low_latitude(1)  mix_low_latitude(1) mix_low_latitude(1)      mix_low_latitude(1)  mix_low_latitude(1) mix_low_latitude(2)   mix_low_latitude(3)  mix_low_latitude(3)   mix_low_latitude(3)   mix_low_latitude(3)  mix_low_latitude(3)      mix_low_latitude(3)    mix_low_latitude(3)]';
  h_l=  [mix_high_latitude(1) mix_high_latitude(1)  mix_high_latitude(1)     mix_high_latitude(1) mix_high_latitude(1) mix_high_latitude(1)  mix_high_latitude(1)     mix_high_latitude(1) mix_high_latitude(2) mix_high_latitude(3) mix_high_latitude(3) mix_high_latitude(3)  mix_high_latitude(3)  mix_high_latitude(3)  mix_high_latitude(3)   mix_high_latitude(2) mix_high_latitude(3)]';
  fra=   [fract(1)              fract(1)              fract(1)                 fract(1)            fract(1)              fract(1)             fract(1)                   fract(1)            fract(1)            fract(2)             fract(3)             fract(3)                fract(3)             fract(3)                fract(3)            fract(3)             fract(3)]';
  tem=   [T_1(1)                 T_1(1)                 T_1(1)                   T_1(1)              T_1(1)                 T_1(3)               T_1(1)                     T_1(1)             T_1(1)           T_1(2)           T_1(3)                    T_1(3)             T_1(1)               T_1(3)                  T_1(3)               T_1(3)             T_1(3)]';
  
  eva=  [evap(1)              evap(1)                evap(1)                evap(3)            evap(1)                evap(1)                evap(1)                 evap(1)               evap(1)                evap(2)        evap(1)                 evap(3)                evap(3)             evap(3)                 evap(3)              evap(3)           evap(3)]';
  riv=[runoff(1)             runoff(2)             runoff(2)             runoff(1)         runoff(1)                runoff(1)              runoff(1)                     runoff(1)        runoff(1)              runoff(2)      runoff(3)             runoff(3)              runoff(3)            runoff(3)                runoff(3)          runoff(3)        runoff(3)]';
    hori=[mix_horiz(1)     mix_horiz(1)               mix_horiz(1)              mix_horiz(1)      mix_horiz(1)              mix_horiz(1)      mix_horiz(2)                   mix_horiz(3)         mix_horiz(1)           mix_horiz(2)     mix_horiz(3)      mix_horiz(3)            mix_horiz(3)         mix_horiz(2)             mix_horiz(1)       mix_horiz(3)     mix_horiz(3)]';
    lohi=[lo_hi_expt(1)   lo_hi_expt(1)             lo_hi_expt(2)              lo_hi_expt(1)    lo_hi_expt(1)               lo_hi_expt(1)    lo_hi_expt(1)                    lo_hi_expt(1)        lo_hi_expt(1)          lo_hi_expt(2)    lo_hi_expt(3)      lo_hi_expt(3)          lo_hi_expt(3)        lo_hi_expt(3)            lo_hi_expt(3)        lo_hi_expt(3)     lo_hi_expt(3)]';
Bering_Strait_expt(ii)=lo_hi_expt(ii);
low=[one(1,1,1,1,1,1,1,1,1,1) one(1,1,1,1,1,1,1,1,1,2)  one(1,1,1,1,1,1,1,1,2,2)   one(1,1,1,1,1,3,1,1,1,1)  one(3,1,1,1,1,1,1,1,1,1) one(1,1,1,3,1,1,1,1,1,1)  one(1,1,1,1,1,2,1,1,1,1) one(1,1,1,1,1,3,1,1,1,1) one(3,1,2,1,1,1,1,1,1,1) one(2,2,2,2,2,2,2,2,2,2) one(3,3,3,3,1,3,3,3,3,3) one(1,3,3,3,3,3,3,3,3,3) one(3,3,3,1,3,3,3,3,3,3) one(3,3,3,3,3,2,3,3,3,3) one(3,3,3,3,3,1,3,3,3,3) one(1,3,2,3,3,3,3,3,3,3) one(3,3,3,3,3,3,3,3,3,3)]';

high=[two(1,1,1,1,1,1,1,1,1,1) two(1,1,1,1,1,1,1,1,1,2)  two(1,1,1,1,1,1,1,1,2,2) two(1,1,1,1,3,1,1,1,1,1) two(3,1,1,1,1,1,1,1,1,1) two(1,1,1,3,1,1,1,1,1,1) two(1,1,1,1,1,2,1,1,1,1) two(1,1,1,1,1,3,1,1,1,1)  two(3,1,2,1,1,1,1,1,1,1) two(2,2,2,2,2,2,2,2,2,2) two(3,3,3,3,1,3,3,3,3,3) two(1,3,3,3,3,3,3,3,3,3) two(3,3,3,1,3,3,3,3,3,3) two(3,3,3,3,3,2,3,3,3,3)  two(3,3,3,3,3,1,3,3,3,3) two(1,3,2,3,3,3,3,3,3,3) two(3,3,3,3,3,3,3,3,3,3) ]';
gr=low-high;
deep=[three(1,1,1,1,1,1,1,1,1,1)  three(1,1,1,1,1,1,1,1,1,2) three(1,1,1,1,1,1,1,1,2,2) three(1,1,1,1,3,1,1,1,1,1) three(3,1,1,1,1,1,1,1,1,1)  three(1,1,1,3,1,1,1,1,1,1) three(1,1,1,1,1,2,1,1,1,1) three(1,1,1,1,1,3,1,1,1,1)  three(3,1,2,1,1,1,1,1,1,1) three(2,2,2,2,2,2,2,2,2) three(3,3,3,3,1,3,3,3,3,3) three(1,3,3,3,3,3,3,3,3,3) three(3,3,3,1,3,3,3,3,3,3) three(3,3,3,3,3,2,3,3,3,3) three(3,3,3,3,3,1,3,3,3,3) three(1,3,2,3,3,3,3,3,3,3) three(3,3,3,3,3,3,3,3,3,3)  ]';
v=[fb(1,1,1,1,1,1,1,1,1,1) fb(1,1,1,1,1,1,1,1,1,2)  fb(1,1,1,1,1,1,1,1,2,2)  fb(1,1,1,1,3,1,1,1,1,1) fb(3,1,1,1,1,1,1,1,1,1) fb(1,1,1,3,1,1,1,1,1,1) fb(1,1,1,1,1,2,1,1,1,1)  fb(1,1,1,1,1,3,1,1,1,1) fb(3,1,2,1,1,1,1,1,1,1) fb(2,2,2,2,2,2,2,2,2,2) fb(3,3,3,3,1,3,3,3,3,3) fb(1,3,3,3,3,3,3,3,3,3) fb(3,3,3,1,3,3,3,3,3,3)  fb(3,3,3,3,3,2,3,3,3,3) fb(3,3,3,3,3,1,3,3,3,3) fb(1,3,2,3,3,3,3,3,3,3)  fb(3,3,3,3,3,3,3,3,3,3) ]';
exps={'all slo','hlf lohi, hlf rr','no lohi, hlf rr','slo: fst evap','slo: fst over','slo: hi-temp','slo: med hor','slo: fast hor','fst-ov md-hl','all med' ,'fst: slo-evap','fst: slo-over','fst: low-temp','fst: med hor' ,'fst: slo hor' ,'fst:slo over, med h_l' ,'all fast'}';
figure (17)
A = table(exps,ov,l_l,h_l,tem,eva,hori,fra, lohi,riv,low,high,deep,v,gr) ;
%z=[4,6,7,8]

table.Border = 'single';
table.ColSep = 'single';
table.RowSep = 'single';
%table.table {border-style: solid; border-bottom-color: rgb(128, 128, 128);border-bottom-width: thin;border-collapse: collapse};


%T(z,:)


% Get the table in string form.
TString = evalc('disp(A)');

% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');

% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');

% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);