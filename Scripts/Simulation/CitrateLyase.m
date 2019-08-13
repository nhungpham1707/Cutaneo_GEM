%% Check ATP citrate liase activity / Identify acetyl-coA source for lipid synthesis 
% % Input: 
% % 
% % - cobra model (model = BuildCryptoGEM); 
% % 
% % Output: 
% % 
% % - Plot with flux for ATP citrate lyase and citrate transporter
% % Simulate:
% % 
% % - Generate new biomass everytime changing Carbon and nitrogen uptake
% %   rate. using func BiomassInDifConditions 
% % - At each C/N ratio constraint the corresponding biomass reaction use in-silico growth rate for C/N < 12.82, and 0.2 for C/N > 12.83 
% % - Change objective function to lipid and optimize for lipid production
% % - Sample the solution space at each C/N ratio
% % - Calculate median for each reaction
% % - Plot results
% % 

%% Nhung 01st May 2019  

%% Nhung 27th June 2019 , modified, used function to create biomass at different C/N 


% ori_model = BuildCryptoGEM;
% 
% clearvars -except ori_model , clc
% model = ori_model; 

clearvars, clc 

load iNP636.mat

%% Set up parameters
gluc = 'r_51_exchange';
urea = 'r_160_exchange';
biomass_N = 'Biomass_nitrogen_abundant';
biomass_Ndel = 'Biomass_nitrogen_deletion'; 
biomass = 'r_1814'; % exchange reaction for biomass
lipid = 'Ex_lipid_body_cytosol';
acyl_N = 'Acyl_Pool_glycerol';
acyl_Ndel = 'Acyl_Pool_glycerol_Ndel';
acyl_glucose = 'Acyl_Pool_glucose'; 
cit_transporter = {'r_1182','r_1183','r_1184', 'r_003_xxx'};

cit_t1 = 'r_1182';
cit_t2 = 'r_1183';
cit_t3 = 'r_1184';
cit_t4 = 'r_003_xxx';
cit_lyase = 'R00352';  % atp-citrate lyase
pdh = 'r_0940'; % pyruvate dehydrogenase
acs = 'r_0127'; % acetyl-coA synthetase
acc = 'r_0124'; % acetyl-coA carboxylase 
mal_acp = 'r_0695'; % malonyl-coa-ACP transcylase 
ptoh = 'R02239_er' ; % phosphatidic acid phosphatase
phos = 'r_25_cl'; % lyso-phosphatidic acid acyltransferase 


model = changeRxnBounds(model,{biomass_Ndel, biomass_N}, 0,'b'); 
model = changeRxnBounds(model,acyl_N,0,'b');
model = changeRxnBounds(model,acyl_Ndel,0,'b');
model = changeRxnBounds(model,{acyl_glucose, acyl_glucose}, [-1000 1000] , {'l','u'});
model = changeRxnBounds(model,{'ATPM'},1,'b');
%% Simulate
CNgg = [0, 10, 12, 14.528, 17.005, 29.882, 32.358, 39.953, 50.024, 71.071];
CNmol = CNgg * 12/14; % 1 g carbon / 1g nitrogen = mol carbon * 12 / mol nitrogen * 14
mol_urea = -25;
model = changeRxnBounds(model,urea,mol_urea,'b');

for i = 1: length(CNmol)
glc_mol (i)  =  CNmol(i)* abs(mol_urea)/6;
end

model = changeObjective(model,biomass);

% get growth rate
for i = 1: length(glc_mol)
    model1 = BiomassInDifConditions(model,-glc_mol(i),6, mol_urea,1);
    model1 = changeRxnBounds(model1,gluc,-glc_mol(i),'b');
    FBA1 = optimizeCbModel(model1);
    gr1 (i) = FBA1.f;
end

%% Sampling 

for i  = 1 : length(CNmol)
    model1 = BiomassInDifConditions(model,-glc_mol(i),6, mol_urea,1); 
    model1 = changeRxnBounds(model1,gluc, - glc_mol(i),'b');
    model_lipid = changeObjective(model1,lipid);

    if CNmol (i) < 12.83
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[gr1(i)*90/100,gr1(i)],{'l','u'});
    else
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[0.2*90/100,0.2],{'l','u'});
    end
    FBA3 = optimizeCbModel(model_lipid);
    lipid3 (i) = FBA3.f; 
    if isempty(FBA3.x) ~= 1
        [sampleStructOut, mixedFraction] = gpSampler (model_lipid,5000,[],120);
        sampling_result{i} = sampleStructOut;
    else
        sampling_result{i} = 'not feasible';
    end
end

%% calculate median flux    
    

check_rxns = {'lipid_synthesis'; acyl_glucose; cit_t1; cit_t2; cit_t3; cit_t4; cit_lyase; pdh; acs; acc; mal_acp; ptoh; phos}; 
check_rxns (:,2) = {'lipid synthesis';'acyl-coA pool';'citrate transport 1' ; 'cit t2';'citrate transport';'cit t4';...
    'ATP-citrate lyase'; 'pyruvate dehydrogenase';'acetyl-CoA synthetase';'acetyl-CoA Carboxylase';...
    'malonyl-CoA-ACP transcylase';'phosphatidic acid phosphatase'; ...
    'lyso-phosphatidic acid acyltransferase'} ; 

plot_rx = {'lipid_synthesis';acyl_glucose;cit_t4; cit_lyase};
plot_rx (:,2) = {'lipid synthesis'; 'acyl-coA pool';'citrate transport from mitochondria to cytosol';'ATP-citrate lyase'}; 
for  i = 1: length(plot_rx)
    index = strmatch(plot_rx(i,1),model.rxns);
    for j = [2:4,6: length(CNmol)] % 0, 5 are infeasible
    median_flux (j,1) = median(sampling_result{j}.points(index,:)); 
    end
    to_plot{i} = median_flux; % 13 column for 13 rxns, 10 value/column for 10 CN ratio
end

% plot 
figure(1)
for i = 1 : length(to_plot)
    plot (CNgg, to_plot{i},'LineWidth',5)
    hold on
end
   ylabel('fluxes mmol/mmol/h')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',30)
    
    legend(plot_rx{1:2,2})
    ylim([0 41])
    xlim([0 71])
    
 figure(3)   
 subplot(1,2,1) 
for i = 1 : length(to_plot)
    plot (CNgg, to_plot{i},'LineWidth',5)
    hold on
      ylabel('fluxes mmol/mmol/h')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',30)
      ylim([0 41])
    xlim([0 71])
    legend(plot_rx{:,2})
end

subplot(1,2,2)
for i = 1: 2
     plot (CNgg, to_plot{i},'LineWidth',5)
    hold on
      ylabel('fluxes mmol/mmol/h')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',30)
%       ylim([0 0.2])
%     xlim([0 71])
    legend(plot_rx{1:2,2})
end
    
    

% legend ('Lipid production','citrate transport 1' , 'cit t2','cit t3', 'cit t4',...
%     'ATP-citrate lyase', 'pyruvate dehydrogenase','acetyl-CoA synthetase','acetyl-CoA Carboxylase',...
%     'malonyl-CoA-ACP transcylase','phosphatidic acid phosphatase', ...
%     'lyso-phosphatidic acid acyltransferase') 
% 
% legend(check_rxns{1:2,2}, check_rxns{6:8,2}, check_rxns{10:12,2})    
    
plotSampleHist(plot_rx(:,1), {sampling_result{4}.points},model,[],[3,7]);

sampleScatterMatrix(plot_rx(:,1),model,sampling_result{4}.points)


%% Plot mean and std with 2 y axis

plot_rx = {'lipid_synthesis';acyl_glucose;cit_t4; cit_lyase};
plot_rx (:,2) = {'lipid synthesis'; 'acyl-coA production';'citrate transport from mitochondria to cytosol';'ATP-citrate lyase'}; 
for  i = 1: length(plot_rx)
    index = strmatch(plot_rx(i,1),model.rxns);
    for j = [2:4,6: length(CNmol)] % 0, 5 are infeasible
    mean_flux (j,1) = mean(sampling_result{j}.points(index,:)); 
    std_flux (j,1) = std(sampling_result{j}.points(index,:));
    end
    mean_plot{i} = mean_flux; % 4 column for 4 rxns, 10 value/column for 10 CN ratio
    std_plot{i} = std_flux;
end


color = {[0 0.4470 0.7410],	[0.8500 0.3250 0.0980], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] }; 
yyaxis left
for i = 1: 2
e = errorbar(CNgg,mean_plot{i},std_plot{i});
% e.Color = color{i};
e.LineWidth =4;
% e.LineStyle = '-';
hold on
end
 ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',40)
       legend(plot_rx(1:2,2))

yyaxis right
for i = 3: 4
e = errorbar(CNgg,mean_plot{i},std_plot{i})
hold on
% e.Color = color{i};
e.LineWidth = 4;
% e.LineStyle = '-';
end
    ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',40)
   xlim([0 72])
legend (plot_rx(:,2))
legend boxoff
box off


figure (8) 
yyaxis left
e = errorbar(CNgg,mean_plot{1},std_plot{1});
e.LineWidth = 4;
e.LineStyle = '--';
 hold on
shadedErrorBar(CNgg,mean_plot{2},std_plot{2},'Lineprops','-b','transparent',1);
% s.mainLine.LineWidth = 6;
    ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',30)
   xlim([0 71])

   yyaxis right
e = errorbar(CNgg,mean_plot{4},std_plot{4});
e.LineWidth = 4;
e.LineStyle = '--';
hold on
% shadedErrorBar(CNgg,mean_plot{3},std_plot{3},'Lineprops','-r','transparent',1);
shadedErrorBar(CNgg,mean_plot{3},std_plot{3},'Lineprops',{'-', 'color',[1 0.49 0] },'transparent',1);
% s2.mainLine.LineWidth = 6;

    ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',30)
   xlim([0 71])
legend(plot_rx{1:4,2})
legend boxoff

box off

figure(3)
yyaxis left
for i = 1:2
    plot(CNgg,mean_plot{i},'LineWidth',5)
    hold on
     ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',30)
   xlim([0 72])
end

yyaxis right 

for i = 3:4
    box off 
    plot(CNgg,mean_plot{i},'LineWidth',5)
    hold on
     ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
    xlabel('C/N (g/g)')
    set(gca,'FontSize',30)
   xlim([0 71])
end

legend(plot_rx{1:4,2})
legend boxoff
box off
% yyaxis left
% %shadedErrorBar(CNgg,mean_plot{1},std_plot{1},'Lineprops','-r','transparent',1);
%  hold on
% shadedErrorBar(CNgg,mean_plot{2},std_plot{2},'Lineprops','-b','transparent',1);
% 
% yyaxis right 
% % shadedErrorBar(CNgg,mean_plot{3},std_plot{3},'Lineprops','-y','transparent',1);
%  hold on
% shadedErrorBar(CNgg,mean_plot{4},std_plot{4},'Lineprops','-g','transparent',1);



%% plot on the same scale with zoom in figure 
plot_rx = {'lipid_synthesis';acyl_glucose;cit_t4; cit_lyase};
plot_rx (:,2) = {'lipid synthesis'; 'acyl-coA production';'citrate transport from mitochondria to cytosol';'ATP-citrate lyase'}; 
for  i = 1: length(plot_rx)
    index = strmatch(plot_rx(i,1),model.rxns);
    for j = [2:4,6: length(CNmol)] % 0, 5 are infeasible
    mean_flux (j,1) = mean(sampling_result{j}.points(index,:)); 
    std_flux (j,1) = std(sampling_result{j}.points(index,:));
    end
    mean_plot{i} = mean_flux; % 4 column for 4 rxns, 10 value/column for 10 CN ratio
    std_plot{i} = std_flux;
end


color = {[0 0.4470 0.7410],	[0.8500 0.3250 0.0980], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] }; 

figure (3) 
subplot(2,1,1)
for i =1: length(plot_rx)
e = errorbar(CNgg,mean_plot{i},std_plot{i});
e.LineWidth = 4;
e.LineStyle = '-';hold on
ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
xlabel('C/N (g/g)')
set(gca,'FontSize',30)
end
xlim([0 72])
max_1 = max(std_plot{3});
Imax_1 = find(std_plot{3} == max_1);
mean_1 = mean_plot{3}; 
ylim([(mean_1(Imax_1) - max_1) (mean_1(Imax_1) + max_1)])
legend (plot_rx(:,2))
legend boxoff
box off
title ('A.')
subplot(2,1,2)
for i = 1:2
e = errorbar(CNgg,mean_plot{i},std_plot{i});
e.LineWidth = 4;
e.LineStyle = '-';hold on
ylabel('Fluxes (mmol.g_{DCW}^{-1}.h^{-1})')
xlabel('C/N (g/g)')
set(gca,'FontSize',30)
end
  xlim([0 72])
  max_2 = max(std_plot{2});
Imax_2 = find(std_plot{2} == max_2);
mean_2 = mean_plot{2}; 
ylim([(mean_2(Imax_2) - max_2) (mean_2(Imax_2) + max_2)])
title ('B.')
box off