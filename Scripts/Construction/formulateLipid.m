function Lipid_model = formulateLipid (draft_model)

%% formulate lipid synthesis in Cryptococcus curvatus. The 5th step in BuildCryptoModel

%% Pham N., 25th May 2018 
model = draft_model;
% Lipid synthesis = 0.9 TAG + 0.07 PE + 0.03 PC . According to Causoni C.
% et al 2017
NewModel = addReaction(model, 'lipid_synthesis','0.9 C00422_er + 0.07 s_1235 + 0.03 s_1230 -> s_1000');
NewModel.rxnNames(end) = cellstr('lipid synthesis');
NewModel.subSystems(strmatch('lipid_synthesis',NewModel.rxns)) = cellstr('Phospholipid synthesis');
% Remove lipid production reaction in yarrowia lypolitica

modelOut = removeRxns(NewModel,NewModel.rxns(strmatch('r_1816',NewModel.rxns))); % remove lipid production reaction
% lipid synthesis pathway is completed! 

Lipid_model = modelOut;