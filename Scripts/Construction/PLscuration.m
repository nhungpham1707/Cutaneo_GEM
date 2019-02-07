function [PL_model, tocheck] = PLscuration(TAG_model)

%% Curate phospholipids synthesis in C.curvatus draft model. PLscuration is the fourth step in BuildCryptoModel function.
%% Nhung, 23rd May 2018

model = TAG_model;
NewModel = model;

%% Modify phospholipid pathway synthesis by removing all wrong reactions in PLs synthesis and adding correct/missing reactions 

%% phospholipid in C.curvatus include PE, PC and PS [Capusoni et al., 2017]. They are synthesised in CDA-DAG and Kennedy pathway. 

%% CDP-DAG pathway 
% CDP-DAG is the main route that produce PE and PC when there is no choline
% and ethanolamine in the growth medium (Henry S.A. et al., 2012) 
% Most of the reactions happen in ER but PS-> PE which happen in mitochondrion (m0) (Carman & Kersting, 2004; Fakas, 2017b). 

% The first step is the synthesis of CDP-DAG from phosphatidate happen in
% ER (Flis V. and Daum G., 2013) -> 'r_0284' 

% The next step is CDP-DAG -> PS in ER (Flis V. and Daum G., 2013) -> 'r_0854'
model.subSystems(strmatch('r_0854',model.rxns)) = cellstr('Phospholipid synthesis'); 

% In the draft model, this step also happen in mitochondria -> 'r_0853', which is not
% consistent with literature. This reaction is removed
% toremove = {'r_0853'}; 
% NewModel = removeRxns(model,toremove);  removal of this reaction lead to
% no growth, leave it here for now. 

% PS in ER is then transported to m0 via a cytosolic factor (Flis V and Daum G., 2013; Junker M
% and Rapoport TA, 2015). It is modeled as PS_ER -> PS_c0 [there is no in the model] ; PS_c0  -> PS_ m0
% ['r_1468' ]
NewModel = addReaction(NewModel,'transport_PS_ER_c0','s_1221 -> s_1220');
NewModel.subSystems(end) = cellstr('Phospholipid synthesis'); 

% PS m0 -> PS c0 and others % check again later
tocheck = ['transport of PS from mito to other compartments'];
% PS -> PE in m0 (Flis V and Daum G,2013) -> 'r_0850'
NewModel.subSystems(strmatch('r_0850',NewModel.rxns)) = cellstr('Phospholipid synthesis'); 
% PE in m0 is transported to ER. This is reversible (Flis V and Daum G., 2013) -> dont have
% that reaction yet. only c0 -> m0 'r_1471'
NewModel = addReaction(NewModel,'transport_PE_m0_ER','s_1237 <=> s_1235');
NewModel.subSystems(end) = cellstr('Phospholipid synthesis'); 


% The next step is a series of methylation of PE to produce PC, happen in
% ER. 
% in the draft model these reactions happen in mitochondrion. This is not
% consistent with literature (Flis V and Daum G, 2013, ect). They are
% removed 
% toremove2 = {'r_0831','r_0874'} ; 
% NewModel = removeRxns(NewModel, toremove2); % leave it here for now

%  in ER, the pathway is not complete
% PE -> PNM 'r_0832' -> PNNM [no reaction] -> PC [no reaction]
NewModel.subSystems(strmatch('r_0832',NewModel.rxns)) = cellstr('Phospholipid synthesis'); 

% S-Adenosyl-L-methionine + Phosphatidyl-N-methylethanolamine <=>
% S-Adenosyl-L-homocysteine + Phosphatidyl-N-dimethylethanolamine 
NewModel = addReaction( NewModel, 'R03424_er','s_1294 + s_1226 -> s_1291 + C04308_er');
NewModel.rxnNames(end) = cellstr('S-adenosyl-L-methionine:phosphatidyl-N-methylethanolamine N-methyltransferase');
NewModel.metNames(strmatch('C04308_er',NewModel.mets)) = cellstr( 'phosphatidyl-N,N-dimethylethanolamine [endoplasmic reticulum]' ); 
NewModel = changeGeneAssociation(NewModel,'R03424_er','YALI0E12441g'); % ssame genes as reactions in mitochondrion r_0831, r_0874  
NewModel.subSystems(end) = cellstr('Phospholipid synthesis'); 
% 	S-Adenosyl-L-methionine + Phosphatidyl-N-dimethylethanolamine <=> S-Adenosyl-L-homocysteine + Phosphatidylcholine
NewModel = addReaction(NewModel,'R01320_er', 's_1294 + C04308_er -> s_1291 + s_1230'); 
NewModel.rxnNames(end)  = cellstr('S-adenosyl-L-methionine:phosphatidyl-N-dimethylethanolamine N-methyltransferase' ); 
NewModel = changeGeneAssociation(NewModel,'R01320_er','YALI0E12441g'); % ssame genes as reactions in mitochondrion r_0831, r_0874  
NewModel.subSystems(end) = cellstr('Phospholipid synthesis'); 
% CDP_DAG pathway is now complete

%% Kennedy pathway - ethanolamine branch 

% The ethanolamine pathway use extracellular ethanolamine from the medium
% ->'r_1249'
NewModel.subSystems(strmatch('r_1249',NewModel.rxns)) = cellstr('Phospholipid synthesis'); 
% check if cryptococcus can synthesize ethanolamine
tocheck = char(tocheck, 'synthesize ethanolamine');

% according to Henry et al, 2012, ethanolamine are phosphorylated with ATP
% in cytosol to produce ethanolamine-P -> 'r_0400'
NewModel.subSystems(strmatch('r_0400',NewModel.rxns)) = cellstr('Phospholipid synthesis'); 
% Ethanolamine-P is further phosphorylated with CTP by enzyme associated with nuclear/ER
% membrane, still at the cytosolic side (Henry et al, 2012) -> 'r_0858'

% The next step CDP-ethanolamine -> PE is happen in ER (Henry et al, 2012). In the draft model,
% this happen in cytosol. It is not clear how CDP-ethanolamine is
% transported to ER. In this model, an artificial transport reaction is
% added, need to check later
tocheck = char(tocheck, 'transport_CDP_ethanolamine_C0_ER');

NewModel = addReaction(NewModel,'transport_CDP_ethanolamine_C0_ER','s_0488 -> C00570_er');
NewModel.subSystems(end) = cellstr('Phospholipid synthesis'); 
NewModel.metNames(end) = cellstr('CDP-ethanolamine [endoplasmic reticulum]');

% 'CDP-ethanolamine [er] + diglyceride [er]  <=> CMP [er] + phosphatidylethanolamine [er] + H+ [er] '
NewModel = addReaction(NewModel, 'R02057_er','C00570_er + s_0597 -> s_0512 + s_1235 + s_0765'); 
NewModel.rxnNames(end) = cellstr('CDPethanolamine:1,2-diacylglycerol ethanolaminephosphotransferase');
NewModel.subSystems(end) = cellstr('Phospholipid synthesis');

% ethanolamine branch is complete!

%% Kennedy pathway - choline branch

% similar to ethanolamine branch, choline is also uptaken from the medium
% -> r_1180
NewModel.subSystems(strmatch('r_1180',NewModel.rxns)) = cellstr('Phospholipid synthesis');

tocheck = char(tocheck, 'synthesize choline');

% choline is also phosphorylated in cytosol to produce choline phosphate by
% ATP. (Henry et al, 2012) -> r_0299
NewModel.subSystems(strmatch('r_0299',NewModel.rxns)) = cellstr('Phospholipid synthesis');
% choline phosphate is phosphorylated to CDP-choline -> r_0300
NewModel.subSystems(strmatch('r_0300',NewModel.rxns)) = cellstr('Phospholipid synthesis');
% CDP-choline need to be transported to ER 
NewModel = addReaction(NewModel,'transport_CDP_choline_C0_ER','s_0483 -> s_0484');
NewModel.subSystems(end) = cellstr('Phospholipid synthesis'); 

% CDP-choline -> PC -> 'r_0303'
NewModel.subSystems(strmatch('r_0303',NewModel.rxns)) = cellstr('Phospholipid synthesis');

%% Choline branch is now complete 

PL_model = NewModel;


