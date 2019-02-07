function Other_model = OtherModification(draft_model)

%% other modification after curate lipid and central carbon metabolism 
NewModel = draft_model;


% Acetyl-CoA carboxylase reaction is not reversible (wikipedia, find other ref) -> r_0124
NewModel.lb (strmatch('r_0124',NewModel.rxns)) = 0;


% ER does not have acetyl-CoA generating enzymes, so it needs to be
% transported from cytosolic by solute carrier family 33, member 1
% (Pietrocola F et al, 2015) 
NewModel = addReaction(NewModel, 'transport_acetyl_coA_C0_ER','s_0380 <=> C00024_er'); 
NewModel.metNames(end) = cellstr('acetyl-CoA [endoplasmic reticulum]' );

% Acetyl-CoA carboxylase is also in ER ref
% https://www.uniprot.org/uniprot/Q00955 
% acetyl-CoA [er] + ATP [er] + bicarbonate [er]  <=> ADP [er] + malonyl-CoA
% [er] + phosphate [er] + H+ [er]
NewModel = addReaction (NewModel, 'R00742_er','C00024_er + s_0447 + C00288_er <=> s_0401 + s_1006 + s_1208 + s_0765'); 
NewModel.rxnNames(end) = cellstr('acetyl-CoA:carbon-dioxide ligase (ADP-forming)');
NewModel.metNames(strmatch('C00288_er',NewModel.mets)) = cellstr('bicarbonate [endoplasmic reticulum]');

% in the draft model, there is no bicarbonate in ER. It is not clear how
% bicarbonate is formed in ER. need to check later
Tocheck = char('formation of bicarbonate');
% 'carbon dioxide [er] + H2O [er]  <=> bicarbonate [er] + H+ [er] '
NewModel = addReaction(NewModel, 'bicarbonate_formation_in_ER','s_0471 + s_1435 <=> C00288_er + s_0765'); 


% add phosphate transport out of c0
NewModel = addReaction (NewModel, 'transport_phosphate_c0_er', ' s_1207 <=> s_1208');

% add transport reaction to move S-adenosyl-L-homocysteine from er to c0 
NewModel = addReaction(NewModel,'transport_S_adenoyl_L_homocysteine_er_c0',{'s_1291','s_1290'},[-1,1],false);
% add transport reaction to move S-adenosyl-L-methionine from c0 to er
NewModel = addReaction(NewModel,'transport_S_adenosyl_L_methionine_c0_er',{'s_1293','s_1294'},[-1,1],false);

Other_model = NewModel;
