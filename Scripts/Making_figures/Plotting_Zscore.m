%% Ploting Zscore from same CN different Cg sampling


TCA_rxns = {'R00352','r_0328', 'r_0330', 'r_0307',  'r_0631',  'r_0629', 'r_1003',  'r_0486', 'r_0488', 'r_0689','r_0804',};
PPP_rxns = {'r_0501', 'r_0106', 'r_0862', 'r_0963','r_0965','r_0880','r_0950','r_1036','r_019_xxx', 'r_0106',  'r_0501', 'r_0040', 'r_0862',  'r_1036', 'r_1035','r_1037','r_0971','r_1035'};
Gly_rxns =  {'r_0573','r_0505', 'r_0105',  'r_0482', 'r_1041', 'r_0865',  'r_0525', 'r_0484', 'r_0859' , 'r_0483', 'r_0866', 'r_0398', 'r_0857', 'r_0941', 'r_0001', 'r_0004', 'r_0181', 'r_0191', 'r_0938', 'r_0127','r_0353'};
lipid = {'Phospholipid synthesis','TAG synthesis','Fatty acid synthesis in cytosol',...
    'Fatty acid desaturation in ER'}

% load reference model from S.cereviae to map amino acids pathway 
cd ('/Users/rubenvanheck/Dropbox/PhD_WUR/Projects/2016_11_20_Cryptococcus/Matlab/Reference_models/S.cerevisiae_model/iMM904 model')
load('iMM904 (2).mat')
ref = iMM904;
ref_formula = printRxnFormula(ref,ref.rxns)%,false, true,true);


load iNP636 
ori_model = model; 



aa = {'S_Alanine_and_Aspartate_Metabolism','S_Arginine_and_Proline_Metabolism',...
'S_Asparagine_metabolism', 'S_Cysteine_Metabolism','S_Glutamate_metabolism',... 
'S_Glutamine_Metabolism', 'S_Glycine_and_Serine_Metabolism', 'S_Histidine_Metabolism',...
'S_Methionine_Metabolism', 'S_Threonine_and_Lysine_Metabolism',...
'S_Tyrosine__Tryptophan__and_Phenylalanine_Metabolism', ...
'S_Valine__Leucine__and_Isoleucine_Metabolism'}; 
    
nu = {'S_Purine_and_Pyrimidine_Biosynthesis'}; 

clc
i = 8;

ref_formula(ismember(ref.subSystems,aa(i)))
    
clc
findComonReaction(model,'pad_c','pac_c')
5-methyltetrahydrofolate(2-)

ala = {'r_0229','r_0233','r_0235','r_0236','r_0237'}; % out of 9 rx in iMM904

arg = {'r_0130','r_0731','r_0217','r_0224','r_0225','r_0226','r_0277','r_0654',...
    'r_0655','r_0656','r_0661','r_0662','r_0659','r_0660','r_0789','r_0790','r_0792','r_0791',...
    'r_0936','r_0657'}; % out of 33
asp = {'r_0229'}; % out of 2 (don't have the one in extracellular space)
cys = {'r_0024' ,'r_0974'}; % out of 10
glu = {'r_0729','r_0296','r_0297','r_0516','r_0508', 'r_0510','r_1004' }; % out of 17
gln = {'r_0221','r_0516','r_0515','r_0513'}; % out of 4
gly = {'r_0174','r_0338','r_0541','r_0546','r_0543','r_0544','r_0545','r_0539',...
    'r_0540', 'r_0538',  'r_0586','r_0587','r_0588','r_0893','r_0664'}; % out of 19
his = {'r_0208','r_0580','r_0575','r_0576','r_0577','r_0604','r_0605','r_0881','r_0882',...
    'r_0008'}; % out of 14
met = { 'r_0160','r_0159','r_0782','r_0783', 'r_0339','r_0340','r_0589','r_0702',...
    'r_0584', 'r_0786'}; % out of 20
threo ={'r_0651','r_0970', 'r_0645','r_0649','r_0537','r_0581','r_0585','r_0033','r_0668',...
    'r_0969','r_0648','r_0667','r_1026','r_1027'}; % out of 19
tyr = {'r_0187','r_0189','r_0203','r_0195','r_0204'};
val = {};

aa_rxns = [ala,arg,cys,glu,gln,his,met,threo,tyr];



% load data Zscore 
cd ('/Users/rubenvanheck/Dropbox/PhD_WUR/Projects/2016_11_20_Cryptococcus/Matlab/2017_11_15_Start_Over/Final_functions/CN ratio analysis/Sampling_increasingC_thesameCN/2019_31_07_sampling_CN6_300_Cmol_162432')

load Zscore_C4.mat
load Zscore_C2.mat
load Zscore_C3.mat
load sampling.mat
Csample{1}.Zscore = Zscore_c2;
Csample{2}.Zscore = Zscore_c3;
Csample{3}.Zscore = Zscore_c4;

CN = [6,8,10, 12,24,30,36,48,60,90,120,180,240, 300] ;
C = [16,24,32] ;



figure(3)

subplot(2,2,1)

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
% legend ('C 16(g)','C 24(g)', 'C 32(g)')
xlabel('C/N (g/g)')
set(gca,'FontSize',30)

% color = {[0.8,0.36,0.36], [0.5,0.5,0.5], [0.67,0.87,0.0],[0 0 0]};
color = {[0.8,0.36,0.36], [0.5,0.5,0.5], [0.2, 0, 1],[0 0 0]};

for i = 1: 3
clear BS x
BS = Csample{i}.Zscore;
x = 1:size(BS,2);     
x_e = x * 2;  
xtick = CN(2:size(BS,2)+1);
subplot(2,2,i+1)
plot(x_e, BS((ismember(model.rxns,TCA_rxns)),:), '.g','MarkerSize',10,...
     'MarkerEdgeColor',color{1},'MarkerFaceColor',color{1});
hold on 
plot(x_e + 0.5,BS((ismember(model.rxns,PPP_rxns)),:), '.r','MarkerSize',10,...
     'MarkerEdgeColor',color{2},'MarkerFaceColor',color{2});
hold on 
plot(x_e +1,BS((ismember(model.rxns,Gly_rxns)),:), '.b','MarkerSize',10,...
     'MarkerEdgeColor',color{3},'MarkerFaceColor',color{3});
hold on 
plot(x_e + 1.5,BS((ismember(model.rxns,aa_rxns)),:), '.','MarkerSize',10,...
    'MarkerEdgeColor',color{4},'MarkerFaceColor',color{4});   
hold on
plot([x,x_e+2],zeros(1,length([x,x_e])),'-k','LineWidth',2)

title (sprintf('%s%4.0f%s','C', C(i),'(g)'))
set(gca, 'XTick',x_e + 1, 'XTickLabel',xtick)             % Set Tick Labels
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     
set(gca,'FontSize',30)
% legend('TCA','PPP','Glycolysis','Amino acids')

end


legend('TCA','PPP','Glycolysis','Amino acids')

%% plot only 3 zscore figure

figure(4)

% color = {[0.8,0.36,0.36], [0.5,0.5,0.5], [0.67,0.87,0.0],[0 0 0]};
color = {[0.8,0.36,0.36], [0.5,0.5,0.5], [0.2, 0, 1],[0 0 0]};
tit = {'A', 'B', 'C'} ;
M_size = 20;
for i = 1: 3
clear BS x
BS = Csample{i}.Zscore;
x = 1:size(BS,2);     
x_e = x * 2;  
xtick = CN(2:size(BS,2)+1);
subplot(3,1,i)
plot(x_e, BS((ismember(model.rxns,TCA_rxns)),:), '.g','MarkerSize',M_size,...
     'MarkerEdgeColor',color{1},'MarkerFaceColor',color{1});
hold on 
plot(x_e + 0.5,BS((ismember(model.rxns,PPP_rxns)),:), '.r','MarkerSize',M_size,...
     'MarkerEdgeColor',color{2},'MarkerFaceColor',color{2});
hold on 
plot(x_e +1,BS((ismember(model.rxns,Gly_rxns)),:), '.b','MarkerSize',M_size,...
     'MarkerEdgeColor',color{3},'MarkerFaceColor',color{3});
hold on 
plot(x_e + 1.5,BS((ismember(model.rxns,aa_rxns)),:), '.','MarkerSize',M_size,...
    'MarkerEdgeColor',color{4},'MarkerFaceColor',color{4});   
hold on
plot([x,x_e+2],zeros(1,length([x,x_e])),'-k','LineWidth',2)

% title (sprintf('%s%4.0f%s','C', C(i),'(g)'))
title (tit (i))
set(gca, 'XTick',x_e + 1, 'XTickLabel',xtick)             % Set Tick Labels
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     
set(gca,'FontSize',30)
% legend('TCA','PPP','Glycolysis','Amino acids')

end

%% Plot each of them on separate plot

M_size = 30; 

i = 3;

figure(i)
BS = Csample{i}.Zscore;
x = 1:size(BS,2);     
x_e = x * 2;  
xtick = CN(2:size(BS,2)+1);
plot(x_e, BS((ismember(model.rxns,TCA_rxns)),:), '.g','MarkerSize',M_size,...
     'MarkerEdgeColor',color{1},'MarkerFaceColor',color{1});
hold on 
plot(x_e + 0.5,BS((ismember(model.rxns,PPP_rxns)),:), '.r','MarkerSize',M_size,...
     'MarkerEdgeColor',color{2},'MarkerFaceColor',color{2});
hold on 
plot(x_e +1,BS((ismember(model.rxns,Gly_rxns)),:), '.b','MarkerSize',M_size,...
     'MarkerEdgeColor',color{3},'MarkerFaceColor',color{3});
hold on 
plot(x_e + 1.5,BS((ismember(model.rxns,aa_rxns)),:), '.','MarkerSize',M_size,...
    'MarkerEdgeColor',color{4},'MarkerFaceColor',color{4});   
hold on
plot([x,x_e+2],zeros(1,length([x,x_e])),'-k','LineWidth',2)

% title (sprintf('%s%4.0f%s','C', C(i),'(g)'))
% title (tit (i))
set(gca, 'XTick',x_e + 1.8, 'XTickLabel',xtick)             % Set Tick Labels
xlabel('C/N (g/g)')                                                 % X Axis Label
ylabel('Zscore')     
set(gca,'FontSize',40)
% legend('TCA','PPP','Glycolysis','Amino acids')

box off

% plot legend
x = [1,3,5,7];
y = [1,1,1,1];
for i = 1: 4
plot( x(i), y(i),'.','MarkerSize',50,...
     'MarkerEdgeColor',color{i},'MarkerFaceColor',color{i})
 hold on
end

set(gca,'visible','off')
text(1.3,1, 'TCA cycle','FontSize', 30)
text(3.3,1,'PPP','FontSize', 30)
text(5.3,1,'Glycolysis','FontSize', 30)
text(7.3,1,'AAs metabolism','FontSize', 30)
xlim([1 8])
ylim([1 1])

%% Checking important reactions 
cit_tca = {'r_0328','r_0330'};
    
formula = printRxnFormula(model,model.rxns,false, true,true);

figure (2)

plot(x, BS((ismember(model.rxns,Gly_rxns)),:))

incre = find(BS((ismember(model.rxns,Gly_rxns)))>0);

dataset(model.rxnNames((ismember(model.rxns,Gly_rxns(incre)))), formula(ismember(model.rxns,Gly_rxns(incre))))

decre = find(BS((ismember(model.rxns,Gly_rxns)))<0);

dataset(model.rxnNames((ismember(model.rxns,Gly_rxns(decre)))), formula(ismember(model.rxns,Gly_rxns(decre))))



test = aa_rxns;

incre = find(BS((ismember(model.rxns,test)))>0);
formula(ismember(model.rxns,test(incre)))
dataset(model.rxnNames((ismember(model.rxns,test(incre)))), formula(ismember(model.rxns,test(incre))))

decre = find(BS((ismember(model.rxns,test)))<0);
formula(ismember(model.rxns,test(decre)))
dataset(model.rxnNames((ismember(model.rxns,test(decre)))), formula(ismember(model.rxns,test(decre))))
