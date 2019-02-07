function C_model = CentralCcuration (Lipid_model)

%% Curate central carbon metabolism in the draft model with correct lipid synthesis from FAcuration, TAGcuration and PLscuration function 
%% CentralCcuration is the last step in curation section in function BuildCryptoModel
%% Nhung, 23rd May 2018

model = Lipid_model;

TCA_rxns = {'r_0328', 'r_0330', 'r_0307',  'r_0631',  'r_0629', 'r_1003',  'r_0486', 'r_0488', 'r_0689','r_0804',};
PPP_rxns = {'r_0501', 'r_0106', 'r_0862', 'r_0963','r_0965','r_0880','r_0950','r_1036','r_019_xxx', 'r_0106',  'r_0501', 'r_0040', 'r_0862',  'r_1036', 'r_1035','r_1037','r_0971','r_1035'};
Gly_rxns =  {'r_0573','r_0505', 'r_0105',  'r_0482', 'r_1041', 'r_0865',  'r_0525', 'r_0484', 'r_0859' , 'r_0483', 'r_0866', 'r_0398', 'r_0857', 'r_0941', 'r_0001', 'r_0004', 'r_0181', 'r_0191', 'r_0938', 'r_0127','r_0353'};
% Glycolysis
model.subSystems(ismember(model.rxns, Gly_rxns)) = cellstr('Glycolysis/gluconeogenesis');    
% Pyruvate to acetyl coA is not in model 1.2.7.11 

% TCA cycle 

% isocitrate to oxalosuccinate is not in here 1.1.1.42
% 2-oxoglutarate (alpha-ketoglutarate) to succinyl-CoA is not in the model
% 1.2.7.11 1.2.7.3
 model.subSystem(ismember(model.rxns,TCA_rxns)) = cellstr('TCA cycle');
 
% PPP
model.subSystem(ismember(model.rxns,PPP_rxns)) = cellstr('Pentose phosphate pathway');

%% change reaction direction for those with wrong directionality 

 
TCA_index = ismember(model.rxns,TCA_rxns);
Gly_index = ismember(model.rxns,Gly_rxns); 
PPP_index = ismember(model.rxns,PPP_rxns);
formula = printRxnFormula(model,model.rxns,false, true,true);

n = 0;
for i = 1: length(TCA_rxns)
    if model.lb(ismember(model.rxns,TCA_rxns(i))) ~= 0 
        n = n+1;
        TCA_reversible (n,1) = TCA_rxns(i);
        TCA_reversible (n,2) = formula(ismember(model.rxns,TCA_rxns(i)));
        TCA_reversible (n,3) = num2cell(model.lb(ismember(model.rxns,TCA_rxns(i))));
    else
    end
end
m = 0;
for i = 1: length(Gly_rxns)
    if model.lb(ismember(model.rxns,Gly_rxns(i))) ~=0 
        m = m+1;
        Gly_rev (m,1) = Gly_rxns(i);
        Gly_rev (m,2) = formula(ismember(model.rxns,Gly_rxns(i)));
        Gly_rev (m,3) = num2cell(model.lb(ismember(model.rxns,Gly_rxns(i))));
    else
    end
end

k = 0;
for i = 1: length(PPP_rxns)
    if model.lb(ismember(model.rxns,PPP_rxns(i))) ~=0
        k = k+1;
        PPP_rev (k,1) = PPP_rxns(i);
        PPP_rev (k,2) = formula(ismember(model.rxns,PPP_rxns(i)));
        PPP_rev (k,3) = num2cell(model.lb(ismember(model.rxns,PPP_rxns(i))));
    else
    end
end

% literature reversible reactions
Lit_rev_TCA = {'r_0330','r_0307','r_0486','r_0629','r_0631','r_0689','r_1003'}'; % 7 reactions 
Lit_rev_Gly = {'r_0398','r_0484','r_0505','r_0525','r_0865','r_0866','r_1041'}'; % 7 reactions
Lit_rev_PPP = {'r_0965','r_1035','r_1036','r_1037'}'; % 4 reactions

% constrast model vs literature
wrong_TCA = setdiff(Lit_rev_TCA,TCA_reversible(:,1))
dataset(model.rxns(ismember(model.rxns,wrong_TCA)),formula(ismember(model.rxns,wrong_TCA)), model.lb(ismember(model.rxns,wrong_TCA)))  % isocitrate -> 2og is not reversible in the model 
% fix TCA 
model.lb(ismember(model.rxns,wrong_TCA)) = -1000;

wrong_Gly = setdiff(Lit_rev_Gly,Gly_rev(:,1)) % correct
wrong_PPP = setdiff(Lit_rev_PPP,PPP_rev(:,1)) % correct 

 reaction = {'r_0865','r_0866','r_0181','r_1247', 'r_1003','r_0488','r_0963' };
 for i = 1:length(reaction)
     index = strmatch(reaction(i), model.rxns);
     model.S(:,index) = model.S(:,index) * -1 ; 
 end
 
 C_model = model; 