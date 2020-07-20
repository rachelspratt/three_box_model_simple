    
%Simple script reads in workspace variables
%from the various scripts: 'arrow_3_box_main_with_loops[sic]..-.m'
%to evaluate the lowest-cost combination
%of input variables:d18O1,d18O2,d18O3,d18Of, observation data
%outputs: 'input matrix', 'best solution', 'best solution location'
%Lisiecki,Spratt 2017


e=load('matlab_expt_1a_110318.mat')
a=load('matlab_expt_1b_110318.mat')
b=load('matlab_expt_2_110418.mat')
%b=load('matlab_ex_3_ch_2.mat')
%b=load('matlab_4_29_18_expt_3_ch_3.mat')
%b=load('matlab_5.mat')
%b=load('matlab_16.mat')


vapor_d18_value=-10.8% calculated value is too low this is a placeholder
low_lat_obs=0.8363;%proxy data LeGrande & Schmidt North Atlantic
hi_lat_obs=-.0402;%volume weighted average
deep_oc_obs=0.0736;
f_obs=vapor_d18_value;
observational_data=[low_lat_obs,hi_lat_obs,deep_oc_obs,f_obs]

cost_box_1=abs((e.box1(:,:,:,:,:,:))-(low_lat_obs));%cost evaluation for each box
cost_box_2=abs((e.box2(:,:,:,:,:,:))-(hi_lat_obs));
cost_box_3=abs((e.box3(:,:,:,:,:,:))-(deep_oc_obs));
cost_fr=abs((e.boxf(:,:,:,:,:,:))-(f_obs));



cost=cost_box_1+cost_box_2+cost_box_3+cost_fr/1000;%Add the cost for each box for each input value

%see the input vectors


best_soln=min(min(min(min(min(min(cost))))));%finds the minimum value in a 6_D array

ind=find(((cost==best_soln)))%numerical location of the element in the array
index_this=e.box1;
sz=size([e.box1]);%size of array

[I1,I2,I3,I4,I5,I6]=ind2sub(sz,ind);%converts element value in a multi dimensional array to an easily readable indexed ouput
best_solution_location=[I1,I2,I3,I4,I5,I6]%linear index location of best solution
compare_to_length=[length(e.ov),length(e.mix_low_latitude),length(e.mix_high_latitude),length(e.alpha),length(e.at_fraction),length(e.mix_horiz)]

%best_numbers=[I1(ov),I2(mix_low_latitude),I3(mix_high_latitude), I4(alpha),I5(at_fraction),I6(mix_horiz)]
best_box1=e.box1(ind)
best_box2=e.box2(ind)
best_box3=e.box3(ind)
best_boxf=e.boxf(ind)
%map_2_input=(ov,mix_low_latitude,mix_high_latitude,alpha,at_fraction,mix_horiz)

figure(1)
subplot (5,3,1)

[X,Y]=meshgrid(e.at_fraction,e.mix_high_latitude);
Z=squeeze(cost(e.I1,e.I2,:,e.I4,:,e.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Hi-Lat Vert Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on


index_this=a.box1;
sz=size([a.box1]);%size of array

[I1,I2,I3,I4,I5,I6]=ind2sub(sz,ind);%converts element value in a multi dimensional array to an easily readable indexed ouput
best_solution_location=[I1,I2,I3,I4,I5,I6]%linear index location of best solution
compare_to_length=[length(a.ov),length(a.mix_low_latitude),length(a.mix_high_latitude),length(a.alpha),length(a.at_fraction),length(a.mix_horiz)]

%best_numbers=[I1(ov),I2(mix_low_latitude),I3(mix_high_latitude), I4(alpha),I5(at_fraction),I6(mix_horiz)]
best_box1=a.box1(ind)
best_box2=a.box2(ind)
best_box3=a.box3(ind)
best_boxf=a.boxf(ind)



subplot(5,3,2)
[X,Y]=meshgrid(a.at_fraction,a.mix_high_latitude);
Z=squeeze(cost(a.I1,a.I2,:,a.I4,:,a.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Hi-Lat Vert Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on



index_this=b.box1;
sz=size([b.box1]);%size of array

[I1,I2,I3,I4,I5,I6]=ind2sub(sz,ind);%converts element value in a multi dimensional array to an easily readable indexed ouput
best_solution_location=[I1,I2,I3,I4,I5,I6]%linear index location of best solution
compare_to_length=[length(b.ov),length(b.mix_low_latitude),length(b.mix_high_latitude),length(b.alpha),length(b.at_fraction),length(b.mix_horiz)]

%best_numbers=[I1(ov),I2(mix_low_latitude),I3(mix_high_latitude), I4(alpha),I5(at_fraction),I6(mix_horiz)]
best_box1=b.box1(ind)
best_box2=b.box2(ind)
best_box3=b.box3(ind)
best_boxf=b.boxf(ind)

subplot(5,3,3)
[X,Y]=meshgrid(b.at_fraction,b.mix_high_latitude);
Z=squeeze(cost(I1,I2,:,I4,:,I6));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Hi-Lat Vert Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on

index_this=e.box1;
sz=size([e.box1]);%size of array

[I1,I2,I3,I4,I5,I6]=ind2sub(sz,ind);%converts element value in a multi dimensional array to an easily readable indexed ouput
best_solution_location=[I1,I2,I3,I4,I5,I6]%linear index location of best solution
compare_to_length=[length(e.ov),length(e.mix_low_latitude),length(e.mix_high_latitude),length(e.alpha),length(e.at_fraction),length(e.mix_horiz)]

%best_numbers=[I1(ov),I2(mix_low_latitude),I3(mix_high_latitude), I4(alpha),I5(at_fraction),I6(mix_horiz)]
best_box1=e.box1(ind)
best_box2=e.box2(ind)
best_box3=e.box3(ind)
best_boxf=e.boxf(ind)
subplot(5,3,4)
[X,Y]=meshgrid(e.at_fraction, e.mix_low_latitude);
Z=squeeze(cost(e.I1,:,e.I3,e.I4,e.I5,:))';
contourf(X,Y,Z,20), colorbar
xlabel('Low-Lat Vert Mix','FontSize', 18)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 18)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on

subplot(5,3,5)
[X,Y]=meshgrid(a.at_fraction, a.mix_low_latitude);
Z=squeeze(cost(a.I1,:,a.I3,a.I4,a.I5,:))';
contourf(X,Y,Z,20), colorbar
xlabel('Low-Lat Vert Mix','FontSize', 18)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 18)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on

index_this=b.box1;
sz=size([b.box1]);%size of array

[I1,I2,I3,I4,I5,I6]=ind2sub(sz,ind);%converts element value in a multi dimensional array to an easily readable indexed ouput
best_solution_location=[I1,I2,I3,I4,I5,I6]%linear index location of best solution
compare_to_length=[length(b.ov),length(b.mix_low_latitude),length(b.mix_high_latitude),length(b.alpha),length(b.at_fraction),length(b.mix_horiz)]

%best_numbers=[I1(ov),I2(mix_low_latitude),I3(mix_high_latitude), I4(alpha),I5(at_fraction),I6(mix_horiz)]
best_box1=b.box1(ind)
best_box2=b.box2(ind)
best_box3=b.box3(ind)
best_boxf=b.boxf(ind)


subplot(5,3,6)
[X,Y]=meshgrid(b.at_fraction, b.mix_low_latitude);
Z=squeeze(cost(b.I1,:,b.I3,b.I4,b.I5,:))';
contourf(X,Y,Z,20), colorbar
xlabel('Low-Lat Vert Mix','FontSize', 18)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 18)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold off

subplot(5,3,7)
[X,Y]=meshgrid(e.ov,e.mix_high_latitude);
Z=squeeze(cost(:,e.I2,:,e.I4,e.I5,e.I6))';
contourf(X,Y,Z,20), colorbar
xlabel('Low-Lat Vert Mix','FontSize', 18)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 18)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)


subplot(5,3,8)
[X,Y]=meshgrid(a.ov,a.mix_high_latitude);
Z=squeeze(cost(:,a.I2,:,a.I4,a.I5,a.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning(Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Atm Flux (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

subplot(5,3,9)
[X,Y]=meshgrid(b.ov,b.mix_high_latitude);
Z=squeeze(cost(:,b.I2,:,b.I4,b.I5,b.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning(Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Atm Flux (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

subplot(5,3,10)
[X,Y]=meshgrid(e.ov,e.alpha);
Z=squeeze(cost(:,e.I2,e.I3,:,e.I5,e.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Atm Flux (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)


subplot(5,3,11)
[X,Y]=meshgrid(a.ov,a.alpha);
Z=squeeze(cost(:,a.I2,a.I3,:,a.I5,a.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Alpha (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

subplot(5,3,12)
[X,Y]=meshgrid(b.ov,b.alpha);
Z=squeeze(cost(:,b.I2,b.I3,:,b.I5,b.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Alpha (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

subplot(5,3,13)
[X,Y]=meshgrid(e.ov,e.at_fraction);
Z=squeeze(cost(:,e.I2,e.I3,e.I4,:,e.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Atm Flux (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on

subplot(5,3,14)
[X,Y]=meshgrid(a.ov,a.at_fraction);
Z=squeeze(cost(:,a.I2,a.I3,a.I4,:,a.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Atm Flux (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on

subplot(5,3,15)
[X,Y]=meshgrid(b.ov,b.at_fraction);
Z=squeeze(cost(:,b.I2,b.I3,b.I4,:,b.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Atm Flux (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold off

figure(10)
[X,Y]=meshgrid(e.alpha,e.mix_low_latitude);
Z=squeeze(cost(e.I1,:,e.I3,:,e.I5,e.I6));
contourf(X,Y,Z,20), colorbar
xlabel('Low Vertical mix (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Alpha (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

%contour plots
figure(12)
[X,Y]=meshgrid(b.alpha,b.mix_high_latitude);
Z=squeeze(cost(I1,I2,:,:,I5,I6));
contourf(X,Y,Z,20), colorbar
xlabel('High Vertical mix (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Alpha (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

figure (8)%test
[X,Y]=meshgrid(e.ov,e.at_fraction);
Z=squeeze(cost(:,I2,I3,I4,:,I6));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Atm Flux (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on


figure(12)
[X,Y]=meshgrid(e.at_fraction,e.alpha);
Z=squeeze(cost(I1,I2,I3,:,:,I6));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Alpha','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)


figure(13)
[X,Y]=meshgrid(e.at_fraction,e.ov);
Z=squeeze(cost(:,I2,I3,I4,:,I6));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Overturning (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

figure(14)
[X,Y]=meshgrid(e.at_fraction,e.mix_low_latitude);
Z=squeeze(cost(I1,:,I3,I4,:,I6));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Low-Lat Vert Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)


figure(15)
[X,Y]=meshgrid(e.at_fraction,e.mix_high_latitude);
Z=squeeze(cost(I1,I2,:,I4,:,I6));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Hi-Lat Vert Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)


figure(16)
[X,Y]=meshgrid(e.at_fraction,e.mix_horiz);
Z=squeeze(cost(I1,I2,I3,I4,:,:));
contourf(X,Y,Z,20), colorbar
xlabel('Atm Flux (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

figure(17)
[X,Y]=meshgrid(e.ov,e.mix_horiz);
Z=squeeze(cost(I1,I2,I3,I4,:,:));
contourf(X,Y,Z,20), colorbar
xlabel('Overturning (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

figure(18)
[X,Y]=meshgrid(e.mix_low_latitude,e.mix_horiz);
Z=squeeze(cost(I1,:,I3,I4,I5,:))';
contourf(X,Y,Z,20), colorbar
xlabel('Low-Lat Vert Mix','FontSize', 18)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 18)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

figure(19)
[X,Y]=meshgrid(e.mix_high_latitude,e.mix_horiz);
Z=squeeze(cost(I1,I2,:,I4,I5,:))';
contourf(X,Y,Z,20), colorbar
xlabel('Hi-Lat Vert Mix (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

figure(20)
[X,Y]=meshgrid(e.alpha,e.mix_horiz);
Z=squeeze(cost(I1,I2,I3,:,I5,:))';
contourf(X,Y,Z,20), colorbar
xlabel('Alpha','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Horiz Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

[X,Y]=meshgrid(e.mix_high_latitude,e.mix_low_latitude);
Z=squeeze(cost(I1,:,:,I4,I5,I6));
contourf(X,Y,Z,20), colorbar
xlabel('High Vert Mix (Sv)','FontSize', 16)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
ylabel('Low Vert Mix (Sv)','FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
hold on



