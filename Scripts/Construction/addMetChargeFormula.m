function model = addMetChargeFormula(model)

% Add met charge and formula for metabolites in iNP636 based on iSB619
% model, Staphylococcus aureus subsp. aureus N315 from BiGG 
% Becker, Scott A., and Bernhard Ø. Palsson. "Genome-scale reconstruction of the metabolic network in Staphylococcus aureus N315: an initial draft to the two-dimensional annotation." BMC microbiology 5.1 (2005): 8.

% Nhung 2018 - metFormula

% load model
load('iSB619.mat')
% load('21082018_final_crypto_model.mat')
load iNP636
w = regexprep(model.metNames,'\[acyl-carrier protein\]','(acyl-carrier protein)','ignorecase');
w = regexprep(w,'\[[a-z]*\s*[a-z]*\]','','ignorecase') % - remove compartment names

w2 = unique(w)


a = intersect(iSB619.metNames,w2)
a = intersect(iSB619.mets,model.mets)


% remove compartment 
a2 = regexprep(a,'_[a-z]$','','ignorecase');
a3 = unique(a2);
a4 = regexprep(model.mets,'_[a-z]$','','ignorecase');
a5 = regexprep(iSB619.mets,'_[a-z]$','','ignorecase');
for i = 1: length(a3)
    clear index index2
    index = strmatch(a3(i),a4,'exact');
    index2 = strmatch(a3(i),a5,'exact');
    model.metCharge(index) = iSB619.metCharge(index2(1));
    model.metFormulas (index) = cellstr(iSB619.metFormulas(index2(1)));
end
