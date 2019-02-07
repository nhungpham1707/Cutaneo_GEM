function TAG_model = TAGcuration(FA_model); 

%% Curate TAG synthesis in the draft model with correct fatty acid synthesis from FAcuration function. TAGcuration is the third step in function BuildCryptoModel
%% Nhung, 23rd May 2018

model = FA_model;

% step 1 _ formation of glycerol-3-phosphate from DHA in glycolysis in
% cytoplasm -> r_0530
% G3P is then transfered to ER -> r_1496
% step 2 _ formation of lyso-phosphatidic acid from G3P + acyl-CoA ->
% r_0534
model.subSystems(strmatch('r_0534',model.rxns)) = cellstr('TAG synthesis');
% step 3_ formation of phosphatidate -> r_25_cl 

model.subSystems(strmatch('r_25_cl',model.rxns)) = cellstr('TAG synthesis');
% step 4_ formation of diacyglycerol -> r_0371 it happens in cytosol not in
% ER -> according to Sorger D and Daum G, 2003. in S.cereviae, there are genes for this
% reaction in ER, m0 but cytosol. For now, just leave the cytosolic
% reaction here and add reaction for ER. 

NewModel = addReaction(model, 'R02239_er', 's_1217 + s_1435 <=> s_0597 + s_1208 + 2 s_0765');
NewModel = changeGeneAssociation(NewModel,'R02239_er','YALI0B14531g');
NewModel.rxnNames(end) = cellstr('1,2-diacyl-sn-glycerol 3-phosphate phosphohydrolase');
NewModel.subSystems(end) = cellstr('TAG synthesis');
% step 5_ formation of TAG -> r_1039 happen in cytosol and using short-chain fa
% instead of acyl-coA. Remove the wrong reaction in cytosol r_1039 
ModelOut = removeRxns(NewModel,'r_1039');
%replace by R02251 diglyceride [er] + acyl-CoA[er]  <=> triglyceride [er] +
%coenzyme A [endoplasmic reticulum]
NewModel = addReaction(ModelOut, 'R02251_er','s_0597 + s_0386 <=> C00422_er + s_0515');
NewModel.metNames(end) = cellstr('triglyceride [endoplasmic reticulum]');
NewModel.rxnNames(end) = cellstr('acyl-CoA:1,2-diacyl-sn-glycerol O-acyltransferase');
NewModel = changeGeneAssociation(NewModel,'R02251_er','YALI0D07986g'); % for EC 2.3.1.20
NewModel.subSystems(end) = cellstr('TAG synthesis');

% acyl coA is generated from all fatty acid with the same amount, need to
% readjust for cryptoccocus
toremove = NewModel.rxns(strmatch('isa acyl-CoA', NewModel.rxnNames));
ModelOut = removeRxns(NewModel,toremove);
% new acyl-CoA production reaction 
% 0.3 C16:0 coA palmitoyl-CoA [endoplasmic reticulum] + 0.15 C18:0 coA stearoyl coA + 0.5 C18:1 coA oleoyl CoA + 0.05 C18:2 coA  octadec-9-ynoyl-CoA (cis) --> acyl CoAs follow iIN800
NewModel = addReaction(ModelOut, 'Acyl_Pool_glucose', '0.3 s_1188 + 0.15 s_1335 + 0.5 s_1148 + 0.05 C02050_er --> s_0386'); % on glucose medium (Beopolous et al, 2011)
% acyl_pool in glycerol from Jan's experiment
NewModel = addReaction(NewModel,'Acyl_Pool_glycerol','0.472 s_1148 + 0.152 s_1188 + 0.197 s_1335 + 0.178 C02050_er  => s_0386'); % C/N 28
NewModel = addReaction(NewModel,'Acyl_Pool_glycerol_Ndel','0.488 s_1148 + 0.197 s_1188 + 0.182 s_1335 + 0.133 C02050_er  => s_0386'); % C/N 2.8

TAG_model = NewModel; 