%% Simulate lipid formation in 2 conditions C/N < = 12.83 mol/mol (growth) and C/N > 12.83 mol/mol (lipid production)
% 
% Input:
% - Model with 2 biomass at normal nitrogen and nitrogen deletion condition. (model = BuildCryptoGEM);
% 
% Set up condition
% 
% - In nitrogen abundant: biomass_N, acyl_pool_glucose
% - In nitrogen deletion: biomass_Ndel, acyl_pool_glucose as in vitro data
% use glucose as carbon source
% 
% Simulate
% 
% - Mode 1. Using biomass in nitrogen abundant
% - Mode 2. Using biomass in nitrogen deletion
% - Mode 3. Using biomass in nitrogen abundant at low C/N (<= 12.83 mol/mol) and biomass at nitrogen deletion at high C/N ratio (> 12.83)


%% Nhung 02nd May 2019 

%% generate model
ori_model = BuildCryptoGEM;

clearvars -except ori_model , clc 

model = ori_model; 


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

%% growth on glucose 
model = changeRxnBounds(model,{acyl_glucose, acyl_glucose}, [0 1000] , {'l','u'});
model = changeRxnBounds(model,acyl_N,0,'b');
model = changeRxnBounds(model,acyl_Ndel,0,'b');
model = changeRxnBounds(model,{'ATPM'},1,'b');

CNgg = [0, 10, 12, 14.528, 17.005, 29.882, 32.358, 39.953, 50.024, 71.071];
CNmol = CNgg * 12/14; % 1 g carbon / 1g nitrogen = mol carbon * 12 / mol nitrogen * 14

% to calculate glucose mol (uptake rate) 
% 1 mol glucose = 6 mol carbon 
% 1 mol Nh4 = 1 mol nitrogen 
mol_urea = -25;
model = changeRxnBounds(model,urea,mol_urea,'b');

for i = 1: length(CNmol)
glc_mol (i)  =  CNmol(i)* abs(mol_urea)/6;
end

%% Mode 1. Simulate the whole thing with biomass nitrogen abundant 
model1 = changeRxnBounds(model,biomass_Ndel, 0,'b'); 
model1 = changeRxnBounds(model1,biomass_N, 1000,'u');
model1 = changeObjective(model1,biomass);

% get growth rate
for i = 1: length(glc_mol)
    model1 = changeRxnBounds(model1,gluc,-glc_mol(i),'b');
    FBA1 = optimizeCbModel(model1);
    gr1 (i) = FBA1.f;
end

%simulate lipid 
model_lipid = changeObjective(model1,lipid);

for i  = 1: length(CNmol)
    model_lipid = changeRxnBounds(model_lipid,gluc,-glc_mol(i),'b');
    if CNmol (i) < 12.83
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[gr1(i)*90/100,gr1(i)],{'l','u'});
    else
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[0.2*90/100,0.2],{'l','u'});
    end
    FBA1 = optimizeCbModel(model_lipid);
    lipid1 (i) = FBA1.f; 
end

%% Mode 2. Simulate using biomass at nitrogen deletion
model2 = changeRxnBounds(model,biomass_Ndel, 1000,'u'); 
model2 = changeRxnBounds(model2,biomass_N, 0,'b');
model2 = changeObjective(model2,biomass);

% get growth rate
for i = 1: length(glc_mol)
    model2 = changeRxnBounds(model2,gluc,-glc_mol(i),'b');
    FBA2 = optimizeCbModel(model2);
    gr2 (i) = FBA2.f;
end

%simulate lipid 
model_lipid = changeObjective(model2,lipid);

for i  = 1: length(CNmol)
        model_lipid = changeRxnBounds(model_lipid,gluc,-glc_mol(i),'b');
    if CNmol (i) < 12.83
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[gr1(i)*90/100,gr1(i)],{'l','u'});
    else
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[0.2*90/100,0.2],{'l','u'});
    end
    FBA2 = optimizeCbModel(model_lipid);
    lipid2 (i) = FBA2.f; 
end


%% Mode 3. Simulate using 2 biomass reaction 
model3 = model; 

for i  = 3: length(CNmol)
    model3 = changeRxnBounds(model3,biomass_Ndel, 0,'b'); 
    model3 = changeRxnBounds(model3,biomass_N, 0,'b');
    model3 = changeRxnBounds(model3,gluc, - glc_mol(i),'b');
    model_lipid = changeObjective(model3,lipid);

    if CNmol (i) < 12.83
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[gr1(i)*90/100,gr1(i)],{'l','u'});
        model_lipid = changeRxnBounds(model_lipid,biomass_N, 1000, 'u');
    else
        model_lipid = changeRxnBounds(model_lipid,{biomass,biomass},[0.2*90/100,0.2],{'l','u'});
        model_lipid = changeRxnBounds(model_lipid,biomass_Ndel,1000,'u');
    end
    FBA3 = optimizeCbModel(model_lipid);
    lipid3 (i) = FBA3.f; 
%     if isempty(FBA3.x) ~= 1
%         [sampleStructOut, mixedFraction] = gpSampler (model_lipid,5000,[],120);
%         sampling_result{i} = sampleStructOut;
%     else
%         sampling_result{i} = 'not feasible';
%     end
end

    
    
    
    
%% plot 3 mode in_vitro against in_silico
figure (4)
plot(CNgg,lipid1,'LineWidth',3) % C/N ratio and lipid
title(sprintf('In-silico lipid body formation %3.2f %s',FBA3.x(strmatch('ATPM',model.rxns)), 'ATP non growth maintenance in \it C.curvatus'))
ylabel('lipid production rate (mol/mol carbon source)')
xlabel('C/N (g/g)')
set(gca,'FontSize',30)

hold on
plot(CNgg,lipid2,'LineWidth',3)

hold on 
plot(CNgg,lipid3,'--','LineWidth',3)


legend ('Mode 1- biomass\_N','Mode 2- biomass\_Ndel', 'Mode 3- biomass\_N for C/N <= 11 (g/g) and biomass\_Ndel for C/N >11 (g/g)')

ylim([0 0.65])
xlim([0 72])

%% Plot lipid in vitro and in-silico
% convert mmol lipid / mmol glucose to g lipid/ g glucose

Mlipid = 15.6; % from Ykema A, et al., 1988 (Jan thesis p48)
Mgluc = 30; % from Jan thesis p.48
 % 1mmol gluc = 0.3 g gluc
 % 1 mmol lipid = 0.156 g lipid
 % 1 mmmol lipid / 1 mmol gluc = 0.156 / 0.3 = 0.52 g/g 
lipid1_g = lipid1/0.52; % to g lipid/ g glucose
lipid2_g = lipid2/0.52;
lipid3_g = lipid3/0.52;

lipid_invitro = [0,0, 0.023, 0.071, 0.052, 0.135, 0.139,0.133, 0.17, 0.162];% experimental data on glucose page 34 fig 3 Jan thesis
%data from Ykema, A., et al 1986

%% Plot results
figure (12)
plot(CNgg,lipid1_g,'LineWidth',3) % C/N ratio and lipid
title(sprintf('In-silico lipid body formation %3.2f %s',FBA3.x(strmatch('ATPM',model.rxns)), 'ATP non growth maintenance in \it C.curvatus'))
ylabel('lipid production rate (g/g carbon source)')
xlabel('C/N (g/g)')
set(gca,'FontSize',30)

hold on
plot(CNgg,lipid2_g,'LineWidth',3)

hold on 
plot(CNgg,lipid3_g,'--','LineWidth',3)

hold on
sz = 100;
scatter(CNgg,lipid_invitro,sz,'MarkerEdgeColor',[0.4660 0.6740 0.1880],...
              'MarkerFaceColor',[0.4660 0.6740 0.1880],...
              'LineWidth',1.5)

legend ('Mode 1- biomass\_N','Mode 2- biomass\_Ndel', 'Mode 3- biomass\_N for C/N <= 11 (g/g) and biomass\_Ndel for C/N >11 (g/g)', 'Experimental data')

ylim([0 0.8])
xlim([0 72])

%% literature figure 3 Jan thesis p48 Ykema, Adrie, et al. "Optimization of lipid production in the oleaginous yeastApiotrichum curvatum in wheypermeate." Applied microbiology and biotechnology 29.2-3 (1988): 211-218.
% ref_result (:,1) = [0, 6,7, 10,20,30,36,40,45,55,58,80]; % C/N ratio g/g
% this is in whey permeate as carbon source
% ref_result (:,2) = [0,0.0031, 0.0031, 0.025,0.075,0.1125,0.125,0.145,0.145,0.15,0.14,0.160]; % lipid g/g carbon
% experimental data on glucose page 34 fig 3 Jan thesis, from Ykema1986
% paper

lipid_invitro = [0,0, 0.023, 0.071, 0.052, 0.135, 0.139,0.133, 0.17, 0.162];% data from Ykema, A., et al 1986

lipid_invitro2 = smooth(lipid_invitro);
hold on
%plot(results1,lipid_invitro2,'-d','color',[1/256 128/256 181/256],'LineWidth',3) % C/N ratio and lipid
sz = 100;
scatter(results1,lipid_invitro2,sz,'MarkerEdgeColor',[0.4660 0.6740 0.1880],...
              'MarkerFaceColor',[0.4660 0.6740 0.1880],...
              'LineWidth',1.5)
legend ('In-silico results','Experimental data')
ylim([0 0.2])




     
    




















