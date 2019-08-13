function [Final_model] = BuildCryptoGEM;

%% This function generate Cryptococcus cuvartus  model using scaffold based approach (Loira et al, 2012). The scaffold model is iNL895 of Y.lipolytica
%% the workflow is as follow 
% - step 1. Draft a model contains reactions from homologous gene list, essential
% reactions, orphan reactions, exchange reactions and biomass reaction
% (func draftModel)
% - Step 2 - 5. Curate model: lipid synthesis 
%     + Fatty acid synthesis (func FAcuration)
%     + TAG synthesis (func TAGcuration)
%     + PLs synthesis (func PLscuration)
% - Step 6-10. Curate model: central carbon metabolism (func CentralCcuration) 
%     + Glycolysis 
%     + TCA
%     + PPP
%     + Others 
% - Step 11. Model at this stage grow really slow, add reactions from non-homologous genes to recover growth
% - Step 12. Update annotation for the new/correct genome
% - Step 13. Add NGAM reaction as reference model does not have 
% - Step 14. Convert metabolite ID
% -Step 15. ATP citrate lyase reaction
% - Step 16. debugging network, add demand reaction 
% - Step 17. set exchange biomass as objective function
% - Step 18. Formulate biomass
% - Step 19. set condition. only glucose now
% - Step 20. Inspecting unbalance reactions  - this was done directly on iNL895.xml prior to construction, save as modified_iNL895.xml
% - Step 21. Convert gene IDs 
% - Step 22. Remove infeasible loop that generate energy from nothing 

%% Data use for this project include
% - iNL895 sbml file which is modified. In the original model, boundary
% metabolites such as hydrogen, coA, Phospho and ect are not read in mat
% file. Leakage reactions were also fixed
% - List of homologous genes from manual annotation with Maarten + Bart annotation
% - '2019_01st_05_to_make_irreversible.xlsx' : list of reactions whose reverse direction involve in loop 
% - '2019_01st_05_to_make_reversible.xlsx': list of reactions whose reverse
% direction help to rescue growth after removing loop

% Nhung, 23rd May 2018
% Nhung, 6th Febuary 2019 - add removing infeasible loop

%% 1. draft model by calling draftModel function
% this step extract reactions from the homologous genes. Orphan, exchange
% reactions, biomass, Essential reactions and important reactions that
% needed for the model to run 

model = draftModel; % this draft model carry flux (FBA.f > 1e-6) 

%% curate model - Lipid synthesis 

% 2. Fatty acid synthesis 

FA_model = FAcuration (model); 

% 3. TAG synthesis 

TAG_model = TAGcuration (FA_model); % acyl-coA pools are set here

% 4. PLs synthesis 

PL_model = PLscuration (TAG_model); 

% 5. Formulate lipid reaction
Lipid_model = formulateLipid (PL_model); 

%% 6. curate model - central carbon metabolism 
C_model = CentralCcuration (Lipid_model); 

%% 7. other curations 

other_model = OtherModification(C_model);

%% 8. modify exchange reactions
Exc_model = swaplbExc(other_model);


%% 9. Remove wrong reactions
Remove_model = RemoveWrongRxns2(Exc_model);  % change this later Nhung

%% 10. Clean model
Clean_model = CleanModel(Remove_model); 

%% Draftting and curating are now completed. 
%% the final model can be used as an input for the ValidateCryptoModel function
model = Clean_model;

%% 11. Model grow really slow. Add essential reactions from non-homologous genes to recover growth
reaction = {'r_0083','r_0090','r_0091','r_0246','r_0471','r_0742','r_0745','r_2_bb'}; 

modelFileName = 'modified_iNL895.xml'; 
dir = fileparts(which(modelFileName)); 
cd (dir)
ref_model = readCbModel(modelFileName); 
ref_formula = printRxnFormula(ref_model,ref_model.rxns);
ref_formula2 = printRxnFormula(ref_model,ref_model.rxns,false,true,true);
% set glycerol medium
medium = {'r_144_exchange',  'r_136_exchange','r_82_exchange', 'r_150_exchange', 'r_134_exchange', 'r_128_exchange', 'r_77_exchange', 'r_162_exchange' }; % 'r_38_exchange' };
C_source = {'r_51_exchange'}; % glycerol
N_source = {'r_160_exchange'}; % Urea 
[er ur] = findExcRxns(model);
model.lb(er) = 0;
model.lb(ismember(model.rxns,medium)) = -1000;
model.lb(ismember(model.rxns, C_source)) = -20;
model.lb(ismember(model.rxns,N_source)) = -1;
FBA = optimizeCbModel(model)

for i = 1: length(reaction)
toadd = ref_formula(ismember(ref_model.rxns,reaction(i)));
NewModel = addReaction(model,reaction{i},toadd{:});
FBA_new = optimizeCbModel(NewModel); 
New_FBA (i,1) = reaction(i);
New_FBA (i,2) = ref_formula2(ismember(ref_model.rxns,reaction(i)));
New_FBA (i,3) = num2cell(FBA_new.f);
end

% added reactions above make model grow faster
%% check if C.curvatus has these genes
genes(:,1) = ref_model.rxns(ismember(ref_model.rxns,reaction));
genes(:,2) = ref_model.grRules(ismember(ref_model.rxns,reaction));

%% add these reactions for further analysis
New_model = model; 
for i = 1: length(reaction)
toadd = ref_formula(ismember(ref_model.rxns,reaction(i)));
New_model = addReaction(New_model,reaction{i},toadd{:});
New_model.rxnNames(end) = ref_model.rxnNames(ismember(ref_model.rxns,reaction(i)));
New_model = changeGeneAssociation(New_model, reaction(i), ref_model.grRules{ismember(ref_model.rxns,reaction(i))});
end

model = New_model;
model_add = addReaction(model,'lipid_body_formation', 's_1000 -> lipid_body_cytosol');
model_add = addExchangeRxn(model_add,model_add.mets(end),0,1000);

model1 = model_add;

%% 12. Update annotation from Bart 

fileName = '13_08_2018_can_remove_from_Maarten_list.xlsx';
dir = fileparts(which(fileName));
cd (dir)

[num,txt,raw] = xlsread(fileName); 

to_remove = raw;

to_remove2 = setdiff(to_remove,{'r_0859','r_0398','r_0941','r_1035'}); % YALI0F15631 = g5226.t1 in Bart annotation
model = removeRxns(model1,to_remove2);

% curate model
model.metNames(strmatch('lipid_body_cytosol',model.metNames)) = cellstr('lipid_body [cytoplasm]');

%% 13. add NGAM reaction (Nhung added on 25th 01 2019 ) , the reference model
% has this reaction with wrong equation (H[extracellular] instead of
% cytoplasm] , already modified in the sbml file r_0249 
% change the name
% model13 = addReaction (model,'ATPM','s_0446 + s_1434 -> s_0400 + s_0764 + s_1207')
% model13.lb(strmatch('ATPM',model13.rxns)) = 1; % set random value for this % in the model, ATPase,cytosolic  
model13 = model; 
model13.rxns(strmatch('r_0249',model13.rxns)) = cellstr('ATPM');

%% 14. Convert metabolite ID
model14 = ConvertID(model13); 

%% 15. ATP citrate lyase reaction
 model15 = CurateATPCitrateLyase(model14); 

 %% 16. debugging network, add demand reaction 
 model_demand = addDemandReaction(model15,'52381_r');
%% 17. set exchange biomass as objective function
model17 = changeObjective(model_demand,'r_1814',1);

%%18. Formulate biomass
model18  = FormulateBiomass(model17); 

%% 19. set condition. only glucose now
model19 = changeRxnBounds(model18,'Acyl_Pool_glycerol',0,'b');
model19 = changeRxnBounds(model19,'Acyl_Pool_glycerol_Ndel',0,'b');
model19 = changeRxnBounds(model19,'Biomass_nitrogen_deletion',0,'b');

%% 20. Inspecting unbalance reactions  - this was done directly on iNL895.xml prior to construction, save as modified_iNL895.xml

%% 21. Convert gene IDs 
model21 = ConvertCryptoGeneID (model19); 


%% 22. Remove infeasible loop that generate energy from thin air 

model22  = RemoveATPLoop(model21);

model22.description = 'Cryptococcus_curvatus_model.xml';

Final_model = model22; 


