function Remove_model = RemoveWrongRxns(draft_model) 

model = draft_model;
%% remove wrong reactions
tol = 1e-6;
%% remove protein synthesis reactions as they are wrong. Each aa can produce protein, that is not correct
toremove = model.rxns(strmatch('protein production',model.rxnNames));
formula = printRxnFormula(model,model.rxns,false, true,true);
formula (strmatch('protein production',model.rxnNames))
assert (length(toremove) == 35) 
modelOut = removeRxns(model,toremove);

model = modelOut;

toremove2 = model.rxns(strmatch('isa dolichol',model.rxnNames));
assert (length(toremove2) == 9 )
modelOut = removeRxns(model, toremove2);
FBA = optimizeCbModel(modelOut); % still run 

%assert (length(model.rxns) - length(modelOut.rxns) == length(toremove2) == 9)
assert (length(model.rxns) - length(modelOut.rxns) == length(toremove2) == 1)
% assert (FBA.f ~=0)
model = modelOut;

semi_model = model; 
toremove_rxns = {'IPS phospholipase C',...
'isa zymosterol ester',...
 'isa fecosterol ester',....
     'isa phosphatidylcholine',...
    'isa ergosterol ester',...
    'palmitoyl transferase for SNarE YAL014C',....
    'palmitoyl transferase for SNarE YAL030W',...
    'palmitoyl transferase for SNarE YDR468C',...
    'palmitoylation of chs3p',... 
    'palmitoylation of YEL013W' } ; 
% toremove_rxns(i) because cryptococcus dont have C24-C26
 n = 0;
toremove_rxns = toremove_rxns'; 
    for i = 1: length(toremove_rxns)
        toremove1 = model.rxns(strmatch(toremove_rxns(i),model.rxnNames)); 
        formula(strmatch(toremove_rxns(i),model.rxnNames))
        modelOut = removeRxns(model,toremove1);
        FBA = optimizeCbModel(modelOut);
        if FBA.f > tol
            model = modelOut;
        else
            n = n+ 1;
            cannot_remove (n) = toremove_rxns(i);
        end
    end  % all 72 reactions were removed with no effects on biomass 
    
    assert (length(model.rxns) + 72 == length(semi_model.rxns))

 Remove_model = model; 