function Exc_model = swaplbExc (draft_model); 
% nutrients are uptaken if their lb is < 0 and in the reaction they are
% reactant mean model.S <0 
model = draft_model;
%% whoever reaction produce metabolite with model.S > 0 change to model.S < 0 
formula = printRxnFormula(model,model.rxns,false, true,true);
[er ur] = findExcRxns(model); 
Exchange(:,1) = model.rxns(er);
Exchange(:,2) = num2cell(model.lb(er));
Exchange(:,3) = num2cell(model.ub(er));
Exchange(:,4) = formula (er);
for j = 1: length(Exchange)
    met = find(model.S(:,strmatch(Exchange(j,1),model.rxns)) ~= 0);  
    if  model.S(met,strmatch(Exchange(j,1), model.rxns)) == 1 ; 
       model.S(met,strmatch(Exchange(j,1), model.rxns)) = -1 ; 
    else
       disp 'no change'
  end
end

% check formula again if they are correc
formula = printRxnFormula(model,model.rxns,false, true,true);
formula(er) % still not correct
% change lb value
store_lb = model.lb(er); 
store_ub = model.ub(er);
% swap ub and lb for exchange reactions from r_1_exchange
index = strmatch('r_1_exchange',Exchange(:,1));
for i = index: length(Exchange)
    indexN = strmatch(Exchange(i,1),model.rxns); 
    model.lb(indexN) = store_ub(i)* -1; 
    model.ub(indexN) = store_lb(i)*-1;
end
formula = printRxnFormula(model,model.rxns,false, true,true);
dataset(formula(er), model.lb(er), model.ub(er)) % now it is correct

Exc_model = model; 