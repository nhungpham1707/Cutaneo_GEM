% prepare list of stoichiometric inconsistent reaction and metabolite  
% Nhung, 7th feb 2019 
model = BuildCryptoGEM; 
ori_model = model;
formula = printRxnFormula(model,model.rxns,false,true,true); 

% remove exchange reactions
[er ur] = findExcRxns(model);

exchange = model.rxns(er); % 188 out of 1546 total rxns

% isolate transport reactions
% transport reactions are those that transport same metabolite to different
% compartmens

n =0;

for i = 1: length(model.rxns)
    reactant = find(model.S(:,i)<0);
    product = find(model.S(:,i)>0);
    reactant_met = regexprep(model.mets(reactant),'_[a-z]+$','');
    produce_met = regexprep(model.mets(product),'_[a-z]+$','');
    same = intersect(reactant_met,produce_met);
    if ~isempty(same)
        n = n+1;
        transport (n,1) = model.rxns(i);
        transport (n,2) = model.rxnNames(i);
        transport (n,3) = formula(i); % 455 transport reactions 
    else
    end
end

% write to file
% fid = fopen('transport_reaction.tsv','w');
% for i = 1: length(transport)
%     fprintf(fid,'%s\t%s\t%s\n',transport{i,1},transport{i,2},transport{i,3})
% end

% normal reactions
internal_rxns = setdiff(model.rxns,[exchange;transport(:,1)]);


% identify inconsistent reactions 


[inform, m, model2] = checkStoichiometricConsistency(model); 

correct_rxns = model2.rxns(model2.SConsistentRxnBool); % only 193 

% need to check  

to_check = setdiff(model.rxns,[correct_rxns;exchange]); % 1165 including transport reaction

% check if there is any reaction in to_check the same in correct_rxns but
% happen in different compartment
compartment (:,1) = {'[cell envelope]','[cytoplasm]','[extracellular]','[mitochondrion]','[nucleus]','[peroxisome]','[endoplasmic reticulum]',...
    '[vacuole]','[Golgi]','[lipid particle]','[vacuolar membrane]'}; % 
correct_formula = regexprep(formula(ismember(model.rxns,correct_rxns)),compartment{:},'');


% find leakage mode - this step takes some time
[Vp, Yp, statp, Vn, Yn, statn] = findMinimalLeakageMode(model2)
[minLeakMetBool, minLeakRxnBool, minSiphonMetBool, minSiphonRxnBool, leakY, siphonY, statp, statn] = findMinimalLeakageModeMet(model2);

% identify leakage model

positive_leak = model2.rxns(any(minLeakRxnBool)); % 435 reactions

negative_leak = model2.rxns(any(minSiphonRxnBool));
n =0 ; 
mode = {} 
for i = 1: length(model2.rxns)
    index = find(minLeakRxnBool(i,:)~=0);
    if ~isempty(index) 
        n = n+1;
        all_index = [i;index']; % one leakage mode
        A {n} = model2.rxns(all_index);
        Asz(n,:) = size(A{n}); % 519 non-unique mode
    else
    end
end

% try again without cofactors
cofactors = {'h2o_x','h2o_v','h2o_n','h2o_m','h2o_l','h2o_g','h2o_e','h2o_r','h2o_c','h_ce','h_c','h_r','h_r','h_e','h_g','h_l','h_m','h_n','h_x','h_v',...
    '18304_c','18304_m','18304_n','13389_c','13389_r','13389_m','13389_n','13389_x','16908_c','16908_r','16908_m','16908_x','18009_c','18009_r','18009_m',...
    '18009_x','18009_e','16474_c','16474_r','16474_m','16474_x','atp_c','atp_r','atp_g','atp_m','atp_n','atp_x','atp_v','datp_c',...,
    'adp_c','adp_r','adp_g','adp_m','adp_n','adp_x','adp_v','dadp_c','dadp_n','amp_c','amp_m','amp_n','amp_x','damp_c','damp_m'}; 

no_factor_model = model2;
index = find(ismember(model2.mets,cofactors));
for i = 1: length(index)
    reaction = find(no_factor_model.S(index(i),:)~=0);
    no_factor_model.S(index(i),reaction) = 0;
end

[inform, m, nof_model2] = checkStoichiometricConsistency(no_factor_model); 
[minLeakMetBool2, minLeakRxnBool2, minSiphonMetBool2, minSiphonRxnBool2, leakY2, siphonY2, statp2, statn2] = findMinimalLeakageModeMet(nof_model2);


inconsit_rxns = nof_model2.rxns(any(minLeakRxnBool2(1:end,:),2));

n = 0;
for i = 1: length(model2.rxns);
    index = find(minLeakRxnBool2(i,:)~=0);
    if ~ isempty (index)
        n = n+1;
        inconsistent (n) = model2.rxns(i);
    else
    end
end


n = 0; 
for i = 1: length(model2.rxns);
    index = find(minLeakRxnBool(i,:)~=0);
    if ~ isempty (index)
        n = n+1;
        inconsistent_factor (n) = model2.rxns(i);
    else
    end
end


% write to file 

% fid = fopen ('2019_01_14th_Positive_leak_reactions.txt','w+') ; % this is
% the file to check 
% for i = 1: length(positive_leak)
%     fprintf(fid,'%s\t%s\t%s\n',positive_leak{i},model.rxnNames{strmatch(positive_leak(i),model.rxns)},formula{strmatch(positive_leak(i),model.rxns)})
% end

