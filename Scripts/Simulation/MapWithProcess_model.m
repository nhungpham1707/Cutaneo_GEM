%% mapping process model prediction to iNP636 model for the effect of C/N ratio on lipid

%% try with the lipid content equation

% y_content = 69.71 - 2.76 * x1 - 297.026 * x2 + 328.79 * x2^2 + 18.97 * x1 * x2 - 29.14 * x1 * x2^2 + 3.13 * x1^2 * x2 - 3.64 * x1^2 * x2^2;
% 
% y_weight = 7.11 + 0.13*x1 - 42.45*x2 + 50.45*x2^2 + 4.70*x1*x2 - 7.58 * x1 *x2^2 + 0.35 *x1^2*x2 - 0.32 * x1^2 * x2^2;
% input data to test the model 
X1 = [8 8 8 16 16 16 24 24 24 32 32 32]; % C Awad et al experimental data
X2 = [0.13 0.26 0.67 0.13 0.26 0.67 0.13 0.26 0.67 0.13 0.26 0.67]; %N  Awad et al experimental data
CN = [60 30 12 120 60 24 180 90 36 240 120 48];



N = [0.1:0.01:0.8];
C = [1.5 : 0.05 : 8.5];
for i = 1: length(C)
    x1 = C(i);
    for j = 1: length(N)
    x2 = N(j);
    lipid_weight (i,j) = 7.11 + 0.13*x1 - 42.45*x2 + 50.45*x2^2 + 4.70*x1*x2 - 7.58 * x1 *x2^2 + 0.35 *x1^2*x2 - 0.32 * x1^2 * x2^2;
    end
end

figure (2)
s = surf(C,N,lipid_weight','LineStyle','none') % it is correct the y equation 
xlabel('C (g.l^{-1})')
ylabel('N (g.l^{-1})')
zlabel('Lipid weight (g.l^{-1})')
set(gca,'FontSize',30)
zlim([0 15])
xlim([1.5 8.5])
ylim([0.1 0.7])
title ('Response surface model (Awad et al 2019)')

x=get(s,'XData');
y=get(s,'YData');
z=get(s,'ZData');
%%Create vectors out of surface's XData and YData
x=x(1,:);
y=y(1,:);
%%Divide the lengths by the number of lines needed
xnumlines = 10; % 10 lines
ynumlines = 10; % 10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);
%%Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k', 'LineWidth', 0.1);
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k', 'LineWidth', 0.1);
end
hold off


%% Now simulate using iNP636

load iNP636
ori_model = model; 
%% simulate lipdi with zero biomass

gluc = 'r_51_exchange';
urea = 'r_160_exchange'; 
lipid = 'Ex_lipid_body_cytosol';
biomass_N = 'Biomass_nitrogen_abundant';
biomass_Ndel = 'Biomass_nitrogen_deletion'; 
biomass = 'r_1814'; % exchange reaction for biomass
acyl_N = 'Acyl_Pool_glycerol';
acyl_Ndel = 'Acyl_Pool_glycerol_Ndel';
acyl_glucose = 'Acyl_Pool_glucose'; 

model = changeRxnBounds(model,{acyl_glucose, acyl_glucose}, [0 1000] , {'l','u'});
model = changeRxnBounds(model,acyl_N,0,'b');
model = changeRxnBounds(model,acyl_Ndel,0,'b');
model = changeRxnBounds(model,{'ATPM'},1,'b');
model = changeRxnBounds(model,{biomass_N, biomass_Ndel},0,'b');


N_mmol = N*1000/14;
C_mol = C/12;
gluc_mmol = C_mol * 1000/ 6; % 6 carbon for glucose 

si_model = model; 
for i = 1: length(gluc_mmol)
    for j = 1: length(N_mmol)     
        si_model2 = si_model; 
        si_model2 = BiomassInDifConditions(si_model2,-gluc_mmol(i),6,-N_mmol(j),1);
        si_model2 = changeRxnBounds(si_model2,gluc,-gluc_mmol(i),'b');
        si_model2 = changeRxnBounds(si_model2,urea,-N_mmol(j),'b');
        si_model2 = changeObjective(si_model2,biomass);
        growth = optimizeCbModel(si_model2,'max');
        gr = growth.f;
        gr_check (i,j) = gr;
        if gr < 0.2
        si_model2 = changeRxnBounds(si_model2,{biomass,biomass},[gr*90/100 gr],{'l','u'});
        else 
        si_model2 = changeRxnBounds(si_model2,{biomass,biomass},[0.2*90/100 0.2],{'l','u'});
        end
        si_model2 = changeObjective(si_model2,lipid);
        FBA_max = optimizeCbModel(si_model2,'max');
        FBA_min = optimizeCbModel(si_model2,'min');
        lipid_max (i,j) = FBA_max.f;
        lipid_min (i,j) = FBA_min.f;
    end
end


lipid_max_g = lipid_max *  475.14 / 1000;
figure (1)
s = surf(C,N,lipid_max_g','LineStyle','none')
xlabel('C (g.l^{-1})')
ylabel('N (g.l^{-1})')
zlabel('lipid production rate (g.g_{DCW}^{-1}.h^{-1})')
set(gca,'FontSize',30)
title('GEM iNP636 model')
xlim([1.5 8.5])
ylim([0.1 0.7])


x=get(s,'XData');
y=get(s,'YData');
z=get(s,'ZData');
%%Create vectors out of surface's XData and YData
x=x(1,:);
y=y(1,:);
%%Divide the lengths by the number of lines needed
xnumlines = 10; % 10 lines
ynumlines = 10; % 10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);
%%Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k', 'LineWidth', 0.1);
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k', 'LineWidth', 0.1);
end
hold off

%% simulate static biomass 
static_model = changeRxnBounds(model,{biomass_N},1000,'u');

for i = 1: length(gluc_mmol)
    for j = 1: length(N_mmol)     
        model3 = static_model; 
        model3 = changeRxnBounds(model3,gluc,-gluc_mmol(i),'b');
        model3 = changeRxnBounds(model3,urea,-N_mmol(j),'b');
        model3 = changeObjective(model3,biomass);
        growth = optimizeCbModel(model3,'max');
        gr = growth.f;
        gr_check (i,j) = gr;
        if gr < 0.2
        model3 = changeRxnBounds(model3,{biomass,biomass},[gr*90/100 gr],{'l','u'});
        else 
        model3 = changeRxnBounds(model3,{biomass,biomass},[0.2*90/100 0.2],{'l','u'});
        end
        model3 = changeObjective(model3,lipid);
        FBA_max = optimizeCbModel(model3,'max');
        FBA_min = optimizeCbModel(model3,'min');
        lipid_max_static (i,j) = FBA_max.f;
        lipid_min (i,j) = FBA_min.f;
    end
end

lipid_max_g_static = lipid_max_static *  475.14 / 1000;
figure (3)
s = surf(C,N,lipid_max_g_static','LineStyle','none')
xlabel('C (g.l^{-1})')
ylabel('N (g.l^{-1})')
zlabel('lipid production rate (g.g_{DCW}^{-1}.h^{-1})')
set(gca,'FontSize',30)
title('GEM iNP636 model - no dynamic biomass')
xlim([1.5 8.5])
ylim([0.1 0.7])
% cmap = cbrewer('seq','GnBu',200);
% colormap(cmap)
% hcb1=colorbar;


x=get(s,'XData');
y=get(s,'YData');
z=get(s,'ZData');
%%Create vectors out of surface's XData and YData
x=x(1,:);
y=y(1,:);
%%Divide the lengths by the number of lines needed
xnumlines = 10; % 10 lines
ynumlines = 10; % 10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);
%%Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k');
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k');
end
hold off

%% Plot everything on the same figure
figure(5)

subplot(3,1,1) % Awad model
surf(C,N,lipid_weight')
xlabel('C (g.l^{-1})')
ylabel('N (g.l^{-1})')
zlabel('Lipid weight (g.l^{-1})')
set(gca,'FontSize',30)
zlim([0 15])
xlim([1.5 8.5])
ylim([0.1 0.7])
title ('Response surface model (Awad et al 2019)')

subplot(3,1,2) % iNP636 with static biomass
surf(C,N,lipid_max_g_static')
xlabel('C (g.l^{-1})')
ylabel('N (g.l^{-1})')
zlabel('lipid production rate (g.g_{DCW}^{-1}.h^{-1})')
set(gca,'FontSize',30)
title('GEM iNP636 model - no dynamic biomass')
xlim([1.5 8.5])
ylim([0.1 0.7])

subplot(3,1,3) %iNP636 with dynamic biomass
surf(C,N,lipid_max_g')
xlabel('C (g.l^{-1})')
ylabel('N (g.l^{-1})')
zlabel('lipid production rate (g.g_{DCW}^{-1}.h^{-1})')
set(gca,'FontSize',30)
title('GEM iNP636 model - with dynamic biomass')
xlim([1.5 8.5])
ylim([0.1 0.7])
