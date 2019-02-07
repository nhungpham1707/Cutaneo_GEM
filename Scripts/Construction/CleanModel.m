function clean_model  = CleanModel(draft_model)

%% clean model by remove unsed genes and  metabolites 
model = draft_model; 

 model = removeTrivialStoichiometry(model);

 %% remove unused genes
unused_gene = model.genes(find(sum(model.rxnGeneMat, 1) == 0)); % there are 14 unsed genes
      
storeGPR = model.grRules;  

model.rxnGeneMat = [];
model.genes = [];
for i = 1 : length(model.rxns)

    model = changeGeneAssociation(model, model.rxns{i}, storeGPR{i});

end

% check again if there is unused gene
model.genes(find(sum(model.rxnGeneMat, 1) == 0)) % empty
% removed 14 unsued genes

%% remove metabolites that not in any reaction

m= 0;
for i = 1: length(model.mets) 
    index = find(model.S(i,:)~=0);
    if isempty(index) == 1
        m = m+1;
        not_use (m) = model.metNames(i);
    else
    end
end

clean_model = model; 