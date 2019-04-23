function modelCrypto_gene = ConvertCryptoGeneID(model); 

%% this function convert yarrow gene ID to crypto gene ID
%% Ortholog genes are gene that share  the same pfam domain, encoding the same enzyme (EC number)
%%  A combination of the protein signatures, EC prediction, BLAST and manual curation was used to find the orthologues
%% Input
% - model after validation step. This is the last step in the model
% reconstruction process
% - data file: list of ortholog genes ('2019_04_19_Yarrow_to_crypto.xlsx'),
% if there are more ortholog crypto genes found, they can be added to this
% file
%% Output
% - model with cryptococcus curvatus gene ID. All yarrow lipolytica genes are removed

%% Nhung 19th April 2019 
fileName = '2019_04_19_Yarrow_to_crypto.xlsx'; 
[num, txt, raw] = xlsread(fileName,1); 

genelist = raw;

for i = 1: length(model.rxns)
model= changeGeneAssociation(model,model.rxns(i),model.grRules{i},genelist(:,1),genelist(:,2));
end

model_newgene = model;


% remove 'YALI0E34155g', there is no ortholog genes % it only participate in one reaction 'r_0402'
%'(YALI0E34155g or g3158.t1 or g6360.t1 or g6988.t1 or g7645.t1 or g921.t1)'
gene = {'YALI0E34155g'};
reaction = find(model_newgene.rxnGeneMat(:,strmatch(gene,model_newgene.genes))~=0);

model_newgene.rxnGeneMat(reaction,strmatch(gene,model_newgene.genes)) = 0 ;
model_newgene = changeGeneAssociation(model_newgene,model_newgene.rxns(reaction),'g3158.t1 or g6360.t1 or g6988.t1 or g7645.t1 or g921.t1');  

unused_gene = model_newgene.genes(find(sum(model_newgene.rxnGeneMat, 1) == 0)); % there are 14 unsed genes
      
storeGPR = model_newgene.grRules;  

model_newgene.rxnGeneMat = [];
model_newgene.genes = [];

for i = 1 : length(model_newgene.rxns)
    model_newgene = changeGeneAssociation(model_newgene, model_newgene.rxns{i}, storeGPR{i});
end


modelCrypto_gene = model_newgene;  % final model have 636 genes