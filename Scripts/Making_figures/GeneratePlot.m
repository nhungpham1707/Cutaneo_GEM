
%% generate plot, figure, and data for poster COBRA conference 

addpath(genpath('/Users/rubenvanheck/Dropbox/PhD_WUR/Projects/2016_11_20_Cryptococcus/Matlab/2017_11_15_Start_Over/Final_functions')); 

% model = BuildCryptoGEM;
clearvars, clc
load iNP636

ori_model = model;
formula = printRxnFormula(model,model.rxns,false, true,true);

%model = ConvertID(model);

compartment (:,1) = {'[cell envelope]','[cytoplasm]','[extracellular]','[mitochondrion]','[nucleus]','[peroxisome]','[endoplasmic reticulum]',...
    '[vacuole]','[Golgi]','[lipid particle]','[vacuolar membrane]'}; 
compartment (:,2) = {'ce','c','e','m','n','x','r','v','g','l','vm'};

comp(:,1) = regexp(model.mets, '_[a-z]+$','match', 'once'); % take compartments
met_comp = regexprep(comp,'_',''); 

% reaction compartment 

n = 0;
for i = 1: length(model.rxns)
    mets = model.mets(find(model.S(:,i)~=0));
    mets2 = regexp(mets,'_[a-z]+$','match','once');
    mets3 = regexprep(mets2,'_','');
    n = n+1;
    reactions (n,1) = model.rxns(i);
    reactions (n,2) = mets3(1); 
end

% orphan reactions
 orphans = findOrphanRxns(model);
 [D, iOr,iAll] = intersect(orphans,reactions(:,1),'stable');
 orphans(:,2) = reactions(iAll,2);
 length(orphans)
 % non-orphan reactions
 non_orphans (:,1) = setdiff(reactions(:,1),orphans);
 [C, iA,iB] = intersect(non_orphans,reactions(:,1),'stable');
 non_orphans(:,2) = reactions(iB,2);
 

label = compartment(:,2)

for i = 1: length(label)
    label_data (i,1)  = nnz(strcmp(met_comp,label(i))); % mets
    label_data (i,2) = nnz(strcmp(orphans(:,2),label(i))); % orphans rxns
    label_data (i,3) = nnz(strcmp(non_orphans(:,2),label(i))); % non_orphan rxns
end

label_data2 = label_data(1:end,2:3);
figure (2)
hb = bar(label_data2)
set(hb(1),'FaceColor',[60/255 180/255 75/255],'EdgeColor',[60/255 180/255 75/255]) 
set(hb(2),'FaceColor',[0 0 128/255],'EdgeColor',[0 0 128/255]) 
% set(hb(3),'FaceColor',

set(gca,'xTickLabel',label, 'fontsize', 50)
title('Distribution of reactions among compartments')
lgd = legend ('# orphans reactions','# non-orphan reactions')
lgd.FontSize = 60;
width = hb.BarWidth;
for i =1 : length(label_data2(:,1))
    row = label_data2(i,:);
    offset = ((width+ 0.5)/length(row))/2;
    x = linspace(i-offset, i+offset, length(row)); 
    text(x,row,num2str(row'),'vert','bottom','horiz','center', 'FontSize',24);
end

% make separate figures
% mets






%% number of transport reactions. Transport reactions are reactions that involve metabolites from different compartments

n = 0;
for i = 1: length(model.rxns)
    mets = model.mets(find(model.S(:,i)~=0));
    mets2 = regexp(mets,'_[a-z]+$','match','once');
    mets3 = regexprep(mets2,'_','');
    mets4 = unique(mets3); 
    if length(mets4) > 1
    n = n+1;
    transport_reaction (n,1) = model.rxns(i); 
    transport_reaction (n,2) = formula(i);
    else
    end
end

%% number of exchange reactions
[er ur] = findExcRxns(model);
exchange_reaction = model.rxns(er);
length(exchange_reaction)

%% number of unique metabolite
unique_mets = unique(regexprep(model.mets,'_[a-z]+$',''));
%% number of unique reactions
internal_reactions = setdiff(model.rxns,[transport_reaction(:,1);exchange_reaction]);

model2 = model;

model2.mets = regexprep(model.mets,'_[a-z]+$','');

formulas_2 = printRxnFormula(model2);
length(unique(formulas_2)) 

figure (3)
% mets
subplot(2,1,1)
number_met = label_data(1:end,1);
hb = bar(number_met,'FaceColor',[112/255 117/255 113/255],'EdgeColor',[112/255 117/255 113/255],'BarWidth',0.3);
set(gca,'xTickLabel',label, 'fontsize', 50)
title('(A)')
lgd = legend ('# Metabolites');
lgd.FontSize = 40;
lgd.Box = 'off';
set(gca,'YTick',[])
width = hb.BarWidth;
for i =1 : length(number_met)
    row = number_met(i);
    x = i;
    text(x,row,num2str(row'),'vert','bottom','horiz','center', 'FontSize',40);
end
set(gca, 'box', 'off')

ylim([0 700])
% reaction
subplot(2,1,2)
hb = bar(label_data2)
set(hb(1),'FaceColor',[60/255 180/255 75/255],'EdgeColor',[60/255 180/255 75/255]) 
set(hb(2),'FaceColor',[0 0 128/255],'EdgeColor',[0 0 128/255]) 
set(gca,'xTickLabel',label, 'fontsize', 50)
title('(B)')
lgd = legend ('# Orphans reactions','# Gene-associated reactions')
lgd.FontSize = 40;
lgd.Box = 'off';
set(gca,'YTick',[])
width = hb.BarWidth;
for i =1 : length(label_data2(:,1))
    row = label_data2(i,:);
    offset = ((width+ 0.5)/length(row))/2;
    x = linspace(i-offset, i+offset, length(row)); 
    text(x,row,num2str(row'),'vert','bottom','horiz','center', 'FontSize',35);
end

set(gca, 'box', 'off')
