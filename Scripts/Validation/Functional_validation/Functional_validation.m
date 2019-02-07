function [data, graph] = endPointValidationCrypto(draft_model)
% - Validate model 
%     + Growth on different carbon sources 
% Input: 
% - data file for carbon and nitrogen sources 
% Ouput:
% - data table includes FBA.f for growth with carbon and nitrogen sources
% Nhung 07th Feb 2019 

model = draft_model; % draft_model is the model from BuildCryptoGEM func
medium = {'r_144_exchange',  'r_136_exchange','r_82_exchange', 'r_150_exchange', 'r_134_exchange', 'r_128_exchange', 'r_77_exchange', 'r_162_exchange' }; % 'r_38_exchange' };

[num, txt, raw] = xlsread('2018_03_07_CandNsources.xlsx',3);
data = raw; 
% constraint pool acyl coa % use the same acyl pool in glucose for every
% carbon sources
model = changeRxnBounds(model,'Acyl_Pool_glucose',1000,'u');
model = changeRxnBounds(model,'Acyl_Pool_glucose',-1000,'l'); 

model = changeRxnBounds(model,'Acyl_Pool_glycerol',0,'b');
model = changeRxnBounds(model,'Acyl_Pool_glycerol_Ndel',0,'b');

[er ur] = findExcRxns(model);
for i = 2: length(data) % remove header
    Newmodel = changeRxnBounds(model, model.rxns(er), 0,'l');
    Newmodel = changeRxnBounds(Newmodel,medium,-1000,'l'); % unlimited mineral and minor nutrients
    Newmodel = changeRxnBounds(Newmodel,medium, 1000, 'u');% unlimited mineral and minor nutrients
    Newmodel = changeRxnBounds(Newmodel, data(i,2), -10,'l'); % nitrogen source
    Newmodel = changeRxnBounds(Newmodel, data(i,4), -10, 'l'); % carbon source
    FBA = optimizeCbModel(Newmodel);
    data(i,5) = num2cell(FBA.f); 
end

data = data;

figure (1)

%% growth slow on glycerol which is not the case in in vitro. Hence, the model needs to be refined in glycerol. 
%% inspect list of non-homologous genes found reactions in ETC. These reactions should be included. 
