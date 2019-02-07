function Draft_model = draftModel

%% Draft cryptococcus model. draftModel is the 1st step in BuildCryptoModel function 
%% Nhung, 23rd May 2018

cd /Users/rubenvanheck/Dropbox/PhD_WUR/Projects/2016_11_20_Cryptococcus/Matlab/2017_11_15_Start_Over/Final_functions/Data

% load scaffold model
model = readCbModel('modified_iNL895_Ruben.xml');
assert (length(model.rxns) == 2002)
assert (length(model.mets) == 1675)
assert (length(model.genes) == 899)
Ref_model=  model;
% Load homologous genes, identify biomass, orphan, and exchange reactions
[num, txt, raw] = xlsread('cryptococcus.inl895.xlsx');
reactions = raw(:,1);
biomass1 = strmatch('biomass production',model.rxnNames);
biomass2 = strmatch('growth',model.rxnNames);
orphans = findOrphanRxns(model);
er = findExcRxns(model);
exchangeRxns = model.rxns(er);
% identify essential reactions
tol = 1e-6;
RxnRatio = singleRxnDeletion(model);
RxnRatio(isnan(RxnRatio)) = 0; % remove nan
EssentialRxns = model.rxns(RxnRatio<tol);
% identify important reactions in non-homologous list 
ImportantRxns = {'r_0892'}; % added on 20thMay2018 r_0892 '3-phospho-serine [cytoplasm] + H2O [cytoplasm]  <=> L-serine [cytoplasm] + phosphate [cytoplasm] '
%reaction = {'r_0083','r_0090','r_0091','r_0246','r_0471','r_0742','r_0745','r_2_bb'}; % from Non-homologous list (in addFromNonHomologous)

% Combine all these reactions
rxnslist2 = [reactions; model.rxns(biomass1); model.rxns(biomass2);orphans;exchangeRxns; EssentialRxns; ImportantRxns]; %reaction']; 
rxnslist = rxnslist2(2:end); %remove the header 'reaction ID'
rxnslist = unique(rxnslist);
% extract submodel from scaffold model
subModel = extractSubNetwork(model,rxnslist);
model = subModel;

assert(length(model.rxns) == length(rxnslist) -13) % -13 reactions named as metaid_r_xxx are not read in the Ref_model 

FBA = optimizeCbModel(model); %% model still run here
assert (FBA.f~= 0)

% finish drafting model. This step generate a draft for further curation
Draft_model = subModel; 
