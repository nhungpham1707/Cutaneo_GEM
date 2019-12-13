function FA_model = FAcuration(draft_model)
%% Curate fatty acid synthesis in the draft model from draftModel function. FAcuration is the second step in function BuildCryptoModel
%% Nhung, 23rd May 2018

model = draft_model;
%% Fatty acid synthesis in cytoplasm 

% Remove wrong reactions 

% In cytoplasm, there are 11 reactions that synthesise fa by
% coorporate malonyl coA to FA instead of FA-ACP, using the same genes as for fa_CoA synthesis. This is not consistent with literature, as there is no free fatty acid in the cell,
% and fatty acid is synthesized indirectly by incoporate malonyl-ACP to FA_ACP.So these reaction are removed
index1  = strmatch('r_0417',model.rxns);
index2 = strmatch('r_0427',model.rxns);
toremove1 = model.rxns(index1:index2);
assert(length(toremove1) == 11) 
% artificial reactions that convert each of C8, c10, c12, ect to fa are
% removed
index3 = strmatch('r_1732',model.rxns);
index4 = strmatch('r_1803',model.rxns);
toremove2 = model.rxns(index3:index4);
assert (length(toremove2) == 72)

% artificial reactions that convert fatty acid to acyl coA, i.e. C8 coA to acyl CoA in ER are removed
index5 = strmatch('r_1643',model.rxns);
index6 = strmatch('r_1720',model.rxns);
toremove6 = model.rxns(index5:index6);
assert (length(toremove6) == 78)
% long chain fa C24-26 are also removed as Cryptococcus does not produce these
% fas
toremove3 = {'r_0437';'r_0438';'r_1261';'r_1262'; 'r_1265';'r_1266';'r_1807';'r_1808'};
% elongation and desaturation in cytosol and mito are also removed because
% this process happen in ER. 
toremove4 = { 'r_0467','r_0727', 'r_0788','r_0810','r_0994','r_0458','r_0460','r_0461','r_0462','r_0463' };

% fatty acid synthesis in cytosol are based on acyl-CoA. This is not
% consistent to literature where it should use acyl-ACP. 
% These reactions are removed and new reactions use ACP are added
toremove5 = {'r_0429','r_0430','r_0464','r_0465','r_0466'};

% remove all above reactions
toremove = [toremove1;toremove2; toremove3; toremove4'; toremove5'; toremove6];
assert (length(toremove) == length(toremove1) + length(toremove2) + length(toremove3) + length(toremove4) + length(toremove5) + length(toremove6))
assert (length(toremove) == 184)

NewModel = removeRxns(model,toremove);

assert(length(NewModel.rxns) == length(model.rxns) - length(toremove));


% Add reactions for fa synthesis in cytosol

% C8:0 lumped
NewModel = addReaction(NewModel, 'R04959_lumped','s_0378 + 3 s_0478 + 6 s_1096 + 9 s_0764 -> 3 s_0470 +  6 s_1091 + C05752_c + 3 s_1434 + 3 s_1495 ');
NewModel.metNames(strmatch('C05752_c',NewModel.mets)) = cellstr('octanoyl-[acyl-carrier protein] [cytoplasm]');
NewModel = changeGeneAssociation(NewModel, 'R04959_lumped','YALI0E23185g and YALI0B15059g and YALI0C11407g and YALI0B19382g'); % same genes for reaction use coA instead of ACP in the Ref_model
NewModel.rxnNames(end) = cellstr('fatty acyl-ACP synthase (n-C8:0ACP),lumped reaction) in cytosol');
NewModel.subSystems(end) = cellstr('Fatty acid synthesis in cytosol');
% C10:0 ACP
NewModel = addReaction(NewModel, 'R04960_R04534_R04535_R04962', 's_0478 + 2 s_1096 + C05752_c + 2 s_0764 -> s_0470 + C05755_c + 2 s_1091 + s_1434 + s_1495');
NewModel.metNames(strmatch('C05755_c',NewModel.mets)) = cellstr('decanoyl-[acyl-carrier protein] [cytoplasm]');
NewModel.rxnNames(end) = cellstr('fatty acyl-ACP synthase (n-C10:0ACP) in cytosol');
NewModel = changeGeneAssociation(NewModel, 'R04960_R04534_R04535_R04962','YALI0E23185g and YALI0B15059g and YALI0C11407g and YALI0B19382g');
NewModel.subSystems(end) = cellstr('Fatty acid synthesis in cytosol');


% C12:0 ACP
NewModel = addReaction(NewModel,'R04963_R04964_R04965_R04725', 's_0478 + 2 s_1096 + C05755_c + 2 s_0764 -> s_0470 + s_0975 + 2 s_1091 + s_1434 + s_1495');
NewModel.rxnNames(end) = cellstr('fatty acyl-ACP synthase (n-C12:0ACP) in cytosol');
NewModel = changeGeneAssociation(NewModel, 'R04963_R04964_R04965_R04725','YALI0E23185g and YALI0B15059g and YALI0C11407g and YALI0B19382g');
NewModel.subSystems(end) = cellstr('Fatty acid synthesis in cytosol');

% C14:0 ACP
NewModel = addReaction(NewModel,'R04726_R04566_R04568_R04967','s_0478 + 2 s_1096 + s_0975 + 2 s_0764 -> s_0470 + s_1042 + 2 s_1091 + s_1434 + s_1495');
NewModel.rxnNames(end) = cellstr('fatty acyl-ACP synthase (n-C14:0ACP) in cytosol');
NewModel = changeGeneAssociation(NewModel, 'R04726_R04566_R04568_R04967','YALI0E23185g and YALI0B15059g and YALI0C11407g and YALI0B19382g');
NewModel.subSystems(end) = cellstr('Fatty acid synthesis in cytosol');

% C16:0 ACP
NewModel = addReaction(NewModel,'R04968_R04543_R04544_R04970','s_0478 + 2 s_1096 + s_1042 + 2 s_0764 -> s_0470 + s_1185 + 2 s_1091 + s_1434 + s_1495');
NewModel.rxnNames(end) = cellstr('fatty acyl-ACP synthase (n-C16:0ACP) in cytosol');
NewModel = changeGeneAssociation(NewModel, 'R04968_R04543_R04544_R04970','YALI0E23185g and YALI0B15059g and YALI0C11407g and YALI0B19382g');
NewModel.subSystems(end) = cellstr('Fatty acid synthesis in cytosol');

% C16:0 ACP --> C16:0 --> C16:0 coA r_0450 and r_0433

% C16:0 CoA is transported from C0 --> ER for further elongation and desaturation r_1447

% C16:0 CoA --> C18:0 CoA in ER 
 NewModel = addReaction(NewModel,'R10825_R10826_R10827_R10828',' s_1006 + 2 s_1097 + s_1188 + 2 s_0765 -> s_0471 + s_0515 + 2 s_1092 + s_1335 + s_1435');
 NewModel.metNames(strmatch('s_1006',NewModel.mets)) = cellstr('malonyl-CoA [endoplasmic reticulum]');
 NewModel.metNames(strmatch('s_1335',NewModel.mets)) = cellstr('stearoyl-CoA [endoplasmic reticulum]');
 NewModel.rxnNames(strmatch('R10825_R10826_R10827_R10828',NewModel.mets)) = cellstr('fatty acyl-ACP synthase (n-C18:0 CoA) in ER');
 NewModel.subSystems(strmatch('R10825_R10826_R10827_R10828',NewModel.mets))  = cellstr('Fatty acid elongation in ER');
 NewModel = changeGeneAssociation(NewModel,'R10825_R10826_R10827_R10828','YALI0E23185g and YALI0B15059g and YALI0C11407g and YALI0B19382g');
 
 % C16:0 CoA --> C16:1 CoA s_0783 in ER 
 NewModel = addReaction(NewModel,'C16coA_desaturase','s_1097 + s_1161 + s_1188 + s_0765 -> s_0783 + s_1092 + 2 s_1435');
 NewModel.metNames(strmatch('s_0783',NewModel.mets)) = cellstr('hexadec-2-enoyl-CoA [endoplasmic reticulum]');
 NewModel.rxnNames (end) = cellstr('palmitoyl-CoA desaturase (n-C16:0CoA -> n-C16:1CoA)');
 NewModel = changeGeneAssociation(NewModel,'FindReactionID','YALI0C05951g');
 NewModel.subSystems(end) = cellstr('Fatty acid desaturation in ER');

 % C18:0 CoA --> C18:1 s_1148 CoA in ER
NewModel = addReaction(NewModel,'R02222_er','s_1097 + s_1161 + s_1335 + s_0765 -> s_1092 + s_1148 + 2 s_1435');
NewModel.rxnNames(end) = cellstr('stearoyl-CoA desaturase (n-C18:0CoA -> n-C18:1CoA)');
NewModel = changeGeneAssociation(NewModel,'R02222_er','YALI0C05951g');
NewModel.subSystems(end) = cellstr('Fatty acid desaturation in ER');


% C18:1 CoA --> C18:2 CoA C01595 (fatty acid not acyl coA), C02050 (acyl coA)  in ER
NewModel = addReaction(NewModel,'R11043_er','s_1097 + s_1161 + s_1148 + s_0765 -> s_1092 + C02050_er + 2 s_1435');
NewModel.metNames(strmatch('C02050_er',NewModel.mets)) = cellstr('octadec-9-ynoyl-CoA [endoplasmic reticulum]');
NewModel.rxnNames(end) = cellstr('Oleoyl-CoA desaturase (n-C18:1CoA -> n-C18:2CoA)');
NewModel = changeGeneAssociation(NewModel,'R11043_er','YALI0C05951g');
NewModel.subSystems(end) = cellstr('Fatty acid desaturation in ER');



% In cytoplasm the hydrolyase steps that remove coA to make fa are missing.
% these reactions need to be added % fatty_acyl-CoA hydrolase

%C16
NewModel = addReaction(NewModel,'R01274_c0',[ NewModel.mets(strmatch('palmitoyl-CoA [cytoplasm]',NewModel.metNames)),NewModel.mets(strmatch('H2O [cytoplasm]',NewModel.metNames)),NewModel.mets(strmatch('coenzyme A [cytoplasm]',NewModel.metNames)),NewModel.mets(strmatch('palmitate [cytoplasm]',NewModel.metNames)) ],[-1 -1 1 1],false, 0, 1000, 0, 'fatty acid synthesis - hydrolyase','YALI0F14135g','',true); 
NewModel.rxnNames(end) = cellstr('fatty acyl-CoA hydrolase');
NewModel.subSystems(strmatch('R01274_c0',NewModel.rxns)) = cellstr('Fatty acid synthesis in cytosol');

% Fatty acid synthesis in cytoplasm completed

% Fatty acid synthesis in mitochondrion is completed but using wrong
% stoichiometric coefficient for hydrogen. 2 instead of 3 hydrogen are required to
% convert Cn-2 -> Cn fatty_ACP. So they need to be changed. 

FA_m0 = {'r_0455','r_0456','r_0457','r_0459'}; 
hydrogen = strmatch('H+ [mitochondrion]', NewModel.metNames); 
NewModel.S(hydrogen,ismember(NewModel.rxns,FA_m0)) = -2; 

% FA synthesis in mitochondria is now completed and correct


FA_model = NewModel; 