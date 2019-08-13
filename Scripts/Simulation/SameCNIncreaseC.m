%% Increase carbon and nitrogen uptake for the same C/N 
%% Sample flux distribution to study if there is any interesting change 
%% to study regulation

% Nhung 24th July 2019

clearvars, clc 

%% Load model 

load iNP636 
ori_model = model; 

%% set parameters
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


%% Define conditions
CN = [0,2,4,6,8,10, 12,24,30,36,48,60,90,120,180,240, 300, 340, 380, 400] ; % g/g CN ratio obtained from Awad et al 2019 
CNmol = CN * 14 /12 ; % mol/mol CN 
% C uptake rate 
C = [8,16,24,32] ; % g C in the meidum obtained from Awad et al 2019 
Cmol = C/12; 
Cmmol = Cmol * 1000/10 ;  % dilute 10 times to get a reasonable simulation data
gluc_mol = Cmol/6; 
gluc_mmol = gluc_mol*1000/10 ; % g

%% Simulate lipid production
for i = 1: length(gluc_mmol)
    for j  = 1: length(CN)
        urea_mmol (i,j) = Cmmol(i)/CN(j); % get urea_mol for the current CN ratio and gluc_mol
        si_model = model; 
        si_model = BiomassInDifConditions(si_model,-gluc_mmol(i),6,-urea_mmol(i,j),1);
        si_model = changeRxnBounds(si_model,gluc,-gluc_mmol(i),'b');
        si_model = changeRxnBounds(si_model,urea,-urea_mmol(i,j),'b');
        
        % constraint biomass / growth
        si_model = changeObjective(si_model,biomass);
        growth = optimizeCbModel(si_model,'max');
        gr = growth.f;
        gr_check (i,j) = gr;
        si_model = changeRxnBounds(si_model,{biomass,biomass},[0.2*90/100 0.2],{'l','u'}); %maximum gr at CN > 11 gg
        
        % simulate lipid 
        si_model = changeObjective(si_model,lipid);
        FBA_max = optimizeCbModel(si_model,'max');
        lipid_max (i,j) = FBA_max.f;
    end
end

%% Plotting lipid production
figure (8)
% plot lipid
% color = {[0 0 0.4]; [0 0 0.9]; [0 0.4 1]; [0 0.6 0.8]};
color = {[0, 0, 1];[0.8500, 0.3250, 0.0980];[0.4940, 0.1840, 0.5560];[0 0.6 0.8]};
subplot(1,2,1)
for i = 1: length(Cmol)
    hold on
plot(CN,lipid_max (i,:),'-x', 'LineWidth',3,...                                                              
            'Color',color{i},...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',3)
leg {i} = sprintf('Cg %6.1f ',C(i));

end
ylabel('Lipid flux (mmol.g_{DCW}^{-1}.h^{-1})')
legend (leg)
xlabel('C/N (g/g)')
set(gca,'FontSize',30)
title('A.')

% plot growth

subplot(1,2,2)
color2 = {[0.8 0.8 0.8];[0.6,0.6,0.6]; [0.4 0.4 0.4]; [0.2 0.2 0.2]};

for i = 1: length(Cmol)
    hold on
plot(CN,gr_check (i,:),'--x', 'LineWidth',3,...
            'Color',color2{i},...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',3)
leg {i} = sprintf('Cg %6.1f ',C(i));
end
legend (leg)
xlabel('C/N (g/g)')
ylabel('Growth rate (1.h^{-1})')
set(gca,'FontSize',30)
title('B.')

%% Sampling 
clear lipid_max gr gr_check CN Cmol Cmmol gluc_mol gluc_mmol urea_mmol 

% define new range of CN and Cmol for sampling based on the plot result
CN = [6,8,10, 12,24,30,36,48,60,90,120,180,240, 300] ;
CNmol = CN * 14 /12 ; % mol/mol CN 
C = [16,24,32] ; % g C in the meidum obtained from Awad et al 2019 
Cmol = C/12; 
Cmmol = Cmol * 1000 /10;  % dilute 10 times to get a reasonable simulation data
gluc_mol = Cmol/6; 
gluc_mmol = gluc_mol*1000 /10 ; % g

for i = 1: length(gluc_mmol)
    n = 0;
    for j  = 1: length(CNmol)
        clear si_model
        urea_mmol (i,j) = Cmmol(i)/CNmol(j); % get urea_mol for the current CN ratio and gluc_mol
        si_model = model; 
        si_model = BiomassInDifConditions(si_model,-gluc_mmol(i),6,-urea_mmol(i,j),1);
        si_model = changeRxnBounds(si_model,gluc,-gluc_mmol(i),'b');
        si_model = changeRxnBounds(si_model,urea,-urea_mmol(i,j),'b');
        
        % constraint biomass / growth
        si_model = changeObjective(si_model,biomass);
        growth = optimizeCbModel(si_model,'max');
        gr = growth.f;
        gr_check (i,j) = gr;
        si_model = changeRxnBounds(si_model,{biomass,biomass},[0.2*90/100 0.2],{'l','u'}); %maximum gr at CN > 11 gg
        
        % simulate lipid 
        si_model = changeObjective(si_model,lipid);
        FBA_max = optimizeCbModel(si_model,'max');
        FBA_min = optimizeCbModel(si_model,'min');
        lipid_max (i,j) = FBA_max.f;
        lipid_min (i,j) = FBA_min.f;
        if isempty (FBA_max.x) ~= 1 
        [sampleStructOut, mixedFraction] = gpSampler (si_model,5000,[],120);
        carbon{i}.sampling_result{j} = sampleStructOut;
        else 
            n = n +1; 
            infeasible_CN (n,1) = Cmol (i);
            infeasible_CN (n,2) = CN(j);
        end        
    end
end

%% Analyze sampling results
Cmol_feasible = Cmol(2:end);
CN_feasible = CN(3:16); % feasible range 4,6,8,10, 12,24,30,36,48,60,90,120,180,240

%% compare reactions at different CN in the same Cmol
 % C1 is infeasible
 % use 2nd Cg and 4th CN  as reference point 
 n = 1; 
for i = 5:length(CNmol)
    Zscore(:,n) = ZscoreSampling(carbon2{4}, carbon2{i});    
    n = n+1;
end

fid = fopen('2019_07_25_Zscore_difCN_Cmol1_4Cmol_10CN.tsv','w+');
for i = 1: length(si_model.rxns)
    fprintf (fid, '%s\t%s\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',si_model.rxns{i},si_model.rxnNames{i},formula{i}...
        ,Zscore(i,1),Zscore(i,2),Zscore(i,3),Zscore(i,4),Zscore(i,5),Zscore(i,6),Zscore(i,7));
end

%% compare reactions at different Cmol in the same CN

%% Plot mean and standard errors for the lipid reactions to choose the CN and C mol to compare

lipid_index = strmatch(lipid,model.rxns,'exact');
y = [];
e =[];
for i = 1:length(Cmol)
      if i == 1; 
        y (i,9:10) = 0;
        e (i,9:10) = 0;
      else
      end
    for j = 1: size(carbon{i}.sampling_result,2)
        
    y (i,j) = mean(carbon{i}.sampling_result{j}.points(index,:));
    e (i,j) = std(carbon{i}.sampling_result{j}.points(index,:)); 
    end
end
x = 1:length(CNmol)
for i = 1: length(Cmol)
    hold on
errorbar((1:length(CNmol)),y(i,:),e(i,:))
lo = y - e;
hi = y + e;

hp = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], 'r');
hold on;
hl = line(x,y);

set(hp, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
set(hl, 'color', 'r', 'marker', 'x');
end


%% Calculate Zscore 

carbon2 = carbon{1}.sampling_result;
carbon3 = carbon{2}.sampling_result;
carbon4 = carbon{3}.sampling_result;

clear urea_sample c_sample
for i = 1: size(carbon4,2)
    urea_sample (i,1) = carbon4{i}.lb(strmatch(urea,model.rxns));
    c_sample (i,1) = carbon4{i}.lb(strmatch(gluc,model.rxns));
    
end

CN_sample = (c_sample*6./urea_sample)*12./14; % 3 cmol gave sampling results for CN6-120, carbon3 give up to 180, carbon4 up to 240


% Calculate Zscore for each Cmol
for i = 1:3
    for j = 2: size(carbon{i}.sampling_result,2)
        Csample{i}.Zscore(:,j-1) = ZscoreSampling(carbon{i}.sampling_result{1}, carbon{i}.sampling_result{j});
    end
end

% save Zscore
Zscore_c2 = Csample{1}.Zscore;
Zscore_c3 = Csample{2}.Zscore;
Zscore_c4 = Csample{3}.Zscore;
save('Zscore_C2.mat','Zscore_c2')
save('Zscore_C3.mat','Zscore_c3')
save('Zscore_C4.mat','Zscore_c4')

%% asign reactions to plot

TCA_rxns = {'R00352','r_0328', 'r_0330', 'r_0307',  'r_0631',  'r_0629', 'r_1003',  'r_0486', 'r_0488', 'r_0689','r_0804',};
PPP_rxns = {'r_0501', 'r_0106', 'r_0862', 'r_0963','r_0965','r_0880','r_0950','r_1036','r_019_xxx', 'r_0106',  'r_0501', 'r_0040', 'r_0862',  'r_1036', 'r_1035','r_1037','r_0971','r_1035'};
Gly_rxns =  {'r_0573','r_0505', 'r_0105',  'r_0482', 'r_1041', 'r_0865',  'r_0525', 'r_0484', 'r_0859' , 'r_0483', 'r_0866', 'r_0398', 'r_0857', 'r_0941', 'r_0001', 'r_0004', 'r_0181', 'r_0191', 'r_0938', 'r_0127','r_0353'};
lipid = {'Phospholipid synthesis','TAG synthesis','Fatty acid synthesis in cytosol',...
    'Fatty acid desaturation in ER'}

to_plot(:,2) = model.subSystems( ~ cellfun(@isempty,model.subSystems))

to_plot(:,1) = model.rxns(~ cellfun(@isempty,model.subSystems))

gly = model.rxns(find(model.S(strmatch('gly_c',model.mets),:)>0));
%ser = find(model.S(strmatch('ser__L_c',model.mets),:)~=0);
%thre = find(model.S(strmatch('thr__L_c',model.mets),:)~=0);
cys = model.rxns(find(model.S(strmatch('cys__L_c',model.mets),:)>0));
lys = model.rxns(find(model.S(strmatch('lys__L_c',model.mets),:)>0));
argi = model.rxns(find(model.S(strmatch('arg__L_c',model.mets),:)>0));
his = model.rxns(find(model.S(strmatch('his__L_c',model.mets),:)>0));

aa = unique([gly;cys;lys;argi;his]);


to_plot = [TCA_rxns,PPP_rxns,Gly_rxns]
% plot for each Cmol 
BS = Zscore_c3;
x = 1:size(BS,2);

figure(1)
for i = 1: 3
    clear BS x
   BS = Csample{i}.Zscore;
x = 1:size(BS,2);     
xtick = CN(2:size(BS,2)+1);
subplot(2,2,i)
plot(x, BS((ismember(model.rxns,TCA_rxns)),:), '.g','MarkerSize',10);  
hold on 
plot(x,BS((ismember(model.rxns,PPP_rxns)),:), 'dr','MarkerSize',5);  
hold on 
plot(x,BS((ismember(model.rxns,Gly_rxns)),:), '+b','MarkerSize',5);  
hold on 
plot(x,BS((ismember(model.rxns,aa)),:), '.k','MarkerSize',10);  
title (sprintf('%s%4.0f%s','C', C(i),'(g)'))
set(gca, 'XTick',x, 'XTickLabel',xtick)             % Set Tick Labels
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     
set(gca,'FontSize',20)
end


legend('TCA','PPP','Glycolysis','Amino acids')


for i = 1: 3 

figure(i) 

BS = Csample{i}.Zscore;
x = 1:size(BS,2);   
xtick = CN(2:size(BS,2)+1);
subplot(2,2,1)
plot(x, BS((ismember(model.rxns,TCA_rxns)),:), '.r','MarkerSize',20);
hold on
plot(x,BS((ismember(model.rxns,PPP_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,Gly_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,aa)),:), '.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
title ('TCA')
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     
set(gca, 'XTick',x, 'XTickLabel',xtick)             % Set Tick Labels

subplot(2,2,2)
plot(x, BS((ismember(model.rxns,PPP_rxns)),:), '.r','MarkerSize',20);
hold on
plot(x,BS((ismember(model.rxns,TCA_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,Gly_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,aa)),:), '.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
title ('PPP')

set(gca, 'XTick',x, 'XTickLabel',xtick)             % Set Tick Labels
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     
subplot(2,2,3)

plot(x, BS((ismember(model.rxns,Gly_rxns)),:), '.r','MarkerSize',20);
hold on
plot(x,BS((ismember(model.rxns,TCA_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,PPP_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,aa)),:), '.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]); 
title ('Glycolysis')
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     
set(gca, 'XTick',x, 'XTickLabel',xtick)             % Set Tick Labels

subplot(2,2,4)

plot(x, BS((ismember(model.rxns,aa)),:), '.r','MarkerSize',20);
hold on
plot(x,BS((ismember(model.rxns,TCA_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,PPP_rxns)),:),'.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]);  
hold on 
plot(x,BS((ismember(model.rxns,Gly_rxns)),:), '.','MarkerSize',20,'MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.5,0.5,0.5]); 

title ('Amino acids')

set(gca, 'XTick',x, 'XTickLabel',xtick)             % Set Tick Labels
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     

end

