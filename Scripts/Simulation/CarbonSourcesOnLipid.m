%% simulate carbon source on lipid formation 
%% Nhung 2019 01st 07

clearvars, clc
%% load model
load iNP636
ori_model = model; 

%% set up parameters
urea = 'r_160_exchange'; 
lipid = 'Ex_lipid_body_cytosol';
biomass_N = 'Biomass_nitrogen_abundant';
biomass_Ndel = 'Biomass_nitrogen_deletion'; 
biomass = 'r_1814'; % exchange reaction for biomass
acyl_N = 'Acyl_Pool_glycerol';
acyl_Ndel = 'Acyl_Pool_glycerol_Ndel';
acyl_glucose = 'Acyl_Pool_glucose'; 

model = removeRxns(model,{acyl_N, acyl_Ndel, biomass_Ndel});
model = changeRxnBounds(model,{'ATPM'},1,'b');
model = changeRxnBounds(model,biomass_N,0,'b');
model = changeRxnBounds(model,{acyl_glucose,acyl_glucose},[0 1000],{'l','u'});

carbon (:,1) = {'r_51_exchange','r_149_exchange', 'r_54_exchange','r_47_exchange','r_59_exchange','r_71_exchange'}; 
carbon (:,2) = {'Glucose','Sucrose','Xylose','Fructose','Ethanol','Glycerol'};  
carbon_number = [6,12,5,6,2,3]; % number of carbon element on the carbon source

%% Simulate lipid and growth for different C source 
% C_mmol = [0.1 : 5 : 300] ;
% urea_mmol = [0.1 : 1 : 30];

N = [0.1:0.01:0.8]; % new data to fit with the process model 
C = [1.5 : 0.05 : 8.5];

% C_mmol = [0.1 : 2 : 100] ; original data before the meeting on 02nd
% august
% urea_mmol = [0.1 : 1 : 10];
urea_mmol = N*1000/14;
C_mmol = C * 1000/12;

for i = 1: length(carbon)
    for j = 1: length(C_mmol)
        for k = 1: length(urea_mmol)
        Csource_mmol = C_mmol(j)/carbon_number(i);    
        model1 = model; 
        model1 = BiomassInDifConditions(model1,-Csource_mmol,carbon_number(i),-urea_mmol(k),1); 
        model1 = changeRxnBounds(model1,carbon(i),-Csource_mmol,'b');
        model1 = changeRxnBounds(model1,urea,-urea_mmol(k),'b');
        FBA = optimizeCbModel(model1);
        gr = FBA.f;
        gr_check(i).flux(j,k) = gr;
        if gr < 0.2
        model1 = changeRxnBounds(model1,{biomass,biomass},[gr*90/100 gr],{'l','u'});
        else 
        model1 = changeRxnBounds(model1,{biomass,biomass},[0.2*90/100 0.2],{'l','u'});
        end
        model1 = changeObjective(model1,lipid);
        FBA_max = optimizeCbModel(model1,'max');
        FBA_min = optimizeCbModel(model1,'min');
        lipid_max(i).flux(j,k) = FBA_max.f;
        lipid_min(i).flux(j,k) = FBA_min.f;
    end
    end
end

%% Plot 
label = {'A','B','C','D','E','F'};
figure(10)
hFig=gcf;
set(hFig, 'Position', [50 50 300 500]);
for i = 1: length(carbon)
    if i ~=5
        n = 10;
    elseif i == 5
        n = 100;
    end
subplot(3,2,i)
contourf(urea_mmol,C_mmol,lipid_max(i).flux,n,'LineStyle','none');
xlabel('N (mmol.g_{DCW}^{-1}.h^{-1})')
% ylabel('Glucose uptake rate (mmol/gDCW/h)')
ylabel('C (mmol.g_{DCW}^{-1}.h^{-1})')
zlabel('Lipid production rate (mmol.g_{DCW}^{-1}.h^{-1})')
set(gca,'FontSize',20)
title(sprintf('%s . %s', label{i}, carbon {i,2}))
%colormap(gca,flipud(gray))
hcb1=colorbar;
end


figure(3)
hFig=gcf;
set(hFig, 'Position', [50 50 300 500]);
for i = 1: length(carbon)
%         if i ~=5
%         n = 10;
%     elseif i == 5
%         n = 100;
%     end
n = 10;
subplot(3,2,i)
contourf(urea_mmol,C_mmol,gr_check(i).flux,n,'LineStyle','none');
xlabel('N (mmol.g_{DCW}^{-1}.h^{-1})')
% ylabel('Glucose uptake rate (mmol/gDCW/h)')
ylabel('C (mmol.g_{DCW}^{-1}.h^{-1})')
zlabel('In-silico growth rate (1.h^{-1})')
set(gca,'FontSize',20)
%title(sprintf('%d. In-silico growth rate (1.h^{-1}) on %s', i, carbon {i,2}))
title(sprintf('%s. %s', label{i}, carbon {i,2}))
cmap = cbrewer('seq','YlGnBu',200);
colormap(cmap)
hcb1=colorbar;
end

xlim([0.1 9.1])
ylim([0.1 98.1])



%% calculate yield 

% load 2019_03rd_07.mat

clear max_C maxValues source
for i = 1: length(carbon) 
    n =0;
    maxValues(i) = max(lipid_max(i).flux(:));
    for j = 1: length(C_mmol)
       if max(lipid_max(i).flux(j,:)) == maxValues(i);
        n = n+1; 
        max_C (n,i) = C_mmol(j);
       else
       end
    end 
end

yield = maxValues* 475.14 ./ max_C * 12


