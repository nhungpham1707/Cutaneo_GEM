%% analyze lipid synthesis and biomass formation at different C/N ratio
%% Nhung 7th Feb 2019

%% generate model
model = BuildCryptoGEM;
ori_model = model;

%% Fix glucose
gluc = 'r_51_exchange';
urea = 'r_160_exchange';

model = changeRxnBounds(model,'r_51_exchange',-10,'b');

%% plot

FBA = optimizeCbModel(model);
max_g = FBA.f;

model = changeRxnBounds(model,'r_1814',FBA.f,'b'); % 0.3376 with -10 glucose

model = changeObjective(model,'Ex_lipid_body_cytosol');

FBA_l = optimizeCbModel(model) % -> no lipid formation

% phenotype phase plane
[growthRates,shadowPrices1,shadowPrices2] = phenotypePhasePlane(model,'r_51_exchange','r_160_exchange');


% product envelope 
model = changeObjective(model,'r_1814');
model = changeRxnBounds(model,'r_1814',0,'l');
model = changeRxnBounds(model,'r_1814',1000,'u');
model = changeRxnBounds(model,'r_51_exchange',-30,'b');


productionEnvelope(model,'N','-k',model.rxns(strmatch('Ex_lipid_body_cytosol',model.rxns)), model.rxns(strmatch('r_1814',model.rxns)))
title('Growth rate and lipid body formation in Crypto')
ylabel('lipid production (mmol/gDCW/h')
xlabel('Growth rate (g/gDCW/h)')

% plot lipid body formation vs carbon concentration with fix maximum
% biomass and fix nitrogen = -0.5;

% -0.5 urea
model = changeRxnBounds(model,'r_1814',max_g,'b'); % 0.3376 with -10 glucose, -1 urea
model = changeRxnBounds(model,'r_160_exchange',-0.5,'b');
model = changeObjective(model,'Ex_lipid_body_cytosol');

for i = 0:50
    model = changeRxnBounds(model,'r_51_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid4 (i+1,1) = i;
    lipid4 (i+1,2) = FBA_l.f;
end

figure(1)
plot(lipid4(:,1),lipid4(:,2),':','LineWidth',8)

% urea -1
model = changeRxnBounds(model,'r_160_exchange',-1,'b');

for i = 0:50
    model = changeRxnBounds(model,'r_51_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid (i+1,1) = i;
    lipid (i+1,2) = FBA_l.f;
end

hold on
plot(lipid(:,1),lipid(:,2),'LineWidth',8)
title('Lipid body formation vs carbon concentration with maximum growth rate in Crypto')
ylabel('Lipid production (mmol/gDCW/h)')
xlabel('Glucose uptake rate (mmol/gDCW/h)')
set(gca,'FontSize',30)

% -5 urea
model = changeRxnBounds(model,'r_160_exchange',-5,'b');
for i = 0:50
    model = changeRxnBounds(model,'r_51_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid2 (i+1,1) = i;
    lipid2 (i+1,2) = FBA_l.f;
end

hold on 
plot(lipid2(:,1),lipid2(:,2),'LineWidth',8)

% -10 urea
model = changeRxnBounds(model,'r_160_exchange',-10,'b');
for i = 0:50
    model = changeRxnBounds(model,'r_51_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid3 (i+1,1) = i;
    lipid3 (i+1,2) = FBA_l.f;
end

hold on 
plot(lipid3(:,1),lipid3(:,2),'LineWidth',8)
legend('Nitrogen (urea) -0.5 (mmol/gDCW/h)','Nitrogen (urea) -1 (mmol/gDCW/h)', 'Nitrogen (urea) -5 (mmol/gDCW/h)', 'Nitrogen (urea) -10 (mmol/gDCW/h)' )

figure (2) 
plot(lipid4(:,1),lipid4(:,2),':','LineWidth',8)
hold on 
plot(lipid3(:,1),lipid3(:,2),'color',[255/255 130/255 1/255],'LineWidth',8)
title('Lipid body formation vs carbon concentration with maximum growth rate in Crypto')
ylabel('Lipid production (mmol/gDCW/h)')
xlabel('Glucose uptake rate (mmol/gDCW/h)')
set(gca,'FontSize',30)
legend('Nitrogen (urea) -0.5 (mmol/gDCW/h)','Nitrogen (urea) -1 , -5, -10 (mmol/gDCW/h)')

% fix carbon with different nitrogen 

model = changeRxnBounds(model,'r_51_exchange',-10,'b');
for i = 0:50
    model = changeRxnBounds(model,'r_160_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid5 (i+1,1) = i;
    lipid5 (i+1,2) = FBA_l.f;
end

figure (2)
plot(lipid5(:,1),lipid5(:,2),'LineWidth',8)
title('Lipid body formation vs nitrogen concentration with maximum growth rate in Crypto')
ylabel('Lipid production (mmol/gDCW/h)')
xlabel('Nitrogen uptake rate (g/gDCW/h)')
set(gca,'FontSize',30)

model = changeRxnBounds(model,'r_51_exchange',-20,'b');
for i = 0:50
    model = changeRxnBounds(model,'r_160_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid6 (i+1,1) = i;
    lipid6 (i+1,2) = FBA_l.f;
end

hold on
plot(lipid6(:,1),lipid6(:,2),'LineWidth',8)

model = changeRxnBounds(model,'r_51_exchange',-30,'b');
for i = 0:50
    model = changeRxnBounds(model,'r_160_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid7 (i+1,1) = i;
    lipid7 (i+1,2) = FBA_l.f;
end

hold on
plot(lipid7(:,1),lipid7(:,2),'LineWidth',8)

model = changeRxnBounds(model,'r_51_exchange',-40,'b');
for i = 0:50
    model = changeRxnBounds(model,'r_160_exchange',-i,'b');
    FBA_l = optimizeCbModel(model);
    lipid8 (i+1,1) = i;
    lipid8 (i+1,2) = FBA_l.f;
end

hold on
plot(lipid8(:,1),lipid8(:,2),'LineWidth',8)
legend('-10 (mmol/gDCW/h) glucose', '-20 (mmol/gDCW/h) glucose', '-30 (mmol/gDCW/h) glucose', '-40 (mmol/gDCW/h) glucose')

%% plot different C/N ratio 

% 16 C/ 1N 28
model = changeRxnBounds(model,'r_51_exchange',-16,'b');
model = changeRxnBounds(model,'r_160_exchange',-1,'b');
model = changeRxnBounds(model,'r_1814', 0, 'l');
model = changeRxnBounds(model, 'r_1814',1000,'u');
model = changeRxnBounds(model,'r_1814',0.3*10/100,'b');
FBA = optimizeCbModel(model) %-> produce lipid with 10% growth 

for i = 1: 90
    model = changeRxnBounds(model,'r_1814',0.3376*i/100,'b');
    FBA = optimizeCbModel(model);
    lipid (i,1) = FBA.f;
    lipid (i,2) = 0.3376*i/100;
end

figure (4)
plot(lipid(:,1),lipid(:,2))
xlabel('lipid production (mmol/gDCW/h)')
ylabel('growth rate (g/gDCW/h)')
% 8 C/ 5 N 2.8 
model = changeRxnBounds(model,'r_51_exchange',-8,'b');
model = changeRxnBounds(model,'r_160_exchange',-5,'b');
for i = 1: 90
    model = changeRxnBounds(model,'r_1814',0.3376*i/100,'b');
    FBA = optimizeCbModel(model);
    lipid2 (i,1) = FBA.f;
    lipid2 (i,2) = 0.3376*i/100;
end
hold on 
plot(lipid2(:,1),lipid2(:,2))
legend('C/N 28', 'C/N 2.8')


figure (5)
plot(lipid(:,2),lipid(:,1))
hold on 
plot(lipid2(:,2),lipid2(:,1))
ylabel('lipid production (mmol/gDCW/h)')
xlabel('growth rate (g/gDCW/h)')
legend('C/N 28', 'C/N 2.8')

% 16C/5N 
model = changeRxnBounds(model,'r_51_exchange',-16,'b');
model = changeRxnBounds(model,'r_160_exchange',-10,'b');
for i = 1: 90
    model = changeRxnBounds(model,'r_1814',0.3376*i/100,'b');
    FBA = optimizeCbModel(model);
    lipid2 (i,1) = FBA.f;
    lipid2 (i,2) = 0.3376*i/100;
end
figure (4)
plot(lipid2(:,2),lipid2(:,1))
legend('C/N 28', 'C/N 2.8')

% convert 16C/1N and 8C/5N to same carbon 
model = changeRxnBounds(model,'r_51_exchange',-8,'b');
model = changeRxnBounds(model,'r_160_exchange',-5,'b');
for i = 1: 90
    model = changeRxnBounds(model,'r_1814',0.3376*i/100,'b');
    FBA = optimizeCbModel(model);
    lipid (i,1) = FBA.f;
    lipid (i,2) = 0.3376*i/100;
end

figure(4)
plot (lipid(:,2),lipid(:,1))
ylabel('lipid production (mmol/gDCW/h)')
xlabel('growth rate (g/gDCW/h)')

model = changeRxnBounds(model,'r_51_exchange',-8,'b');
model = changeRxnBounds(model,'r_160_exchange',-8/16,'b');
for i = 1: 90
    model = changeRxnBounds(model,'r_1814',0.3376*i/100,'b');
    FBA = optimizeCbModel(model);
    lipid2 (i,1) = FBA.f;
    lipid2 (i,2) = 0.3376*i/100;
end
hold on 
plot (lipid2(:,2),lipid2(:,1))
legend('C/N 2.8', 'C/N 28')



%% 2nd simulation approach, dont fix biomass. to prove that lipid only form when excess carbon and limit nitrogen

model = changeRxnBounds(model,'r_1814',max_g*10/100,'l'); % constraint lb to 10% of maximum growth at -10 glucose, ub is unconstraint
model = changeRxnBounds(model,'r_1814',1000,'u');

model = changeObjective(model,'Ex_lipid_body_cytosol'); 

% change carbon concentration and fix nitrogen 
 [lipids,shadowPrices1,shadowPrices2] = phenotypePhasePlane(model,gluc,urea);
    nPts = 50;
    range1 = 20;
    range2 = 20;
ind1 = linspace(0,range1,nPts);
ind2 = linspace(0,range2,nPts);

figure(4);
surfl(ind1,ind2,lipids);
xlabel(' Glucose uptake rate (mmol/g DW-hr)', 'FontSize',30), ylabel(' Urea uptake rate (mmol/g DW-hr)', 'FontSize',30), zlabel('Lipid body formation rate (mmol/g DW-hr)', 'FontSize',30);
colormap (jet)
