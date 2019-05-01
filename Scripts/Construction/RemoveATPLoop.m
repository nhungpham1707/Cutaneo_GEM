function [model] = RemoveATPLoop(model);

% % This function remove thermodynamic infeasible loop that generate energy from nothing. 
% % This is the step in BuildCryptoGEM function 
% % Input:
% % - model: in mat format
% % - '2019_01st_05_to_make_irreversible.xlsx' : list of reactions whose reverse direction involve in loop 
% % - '2019_01st_05_to_make_reversible.xlsx': list of reactions whose reverse direction help to rescue growth
% % 
% % Output:
% % - model that cannot generate ATP when there is no nutrient 
% % 
% % Nhung, 2019_01st_05 

atp_model = changeObjective(model,'ATPM');
biomass_model = model; 

FBA = optimizeCbModel(biomass_model) % growth with nutrient
assert(FBA.f>1e-6)
FBA = optimizeCbModel(atp_model) % 1000 ATP production with nutrient
assert(FBA.f==1000)
% no nutrient model

[er ur] = findExcRxns(model);
medium = model.rxns(ur);
medium_bound = model.lb(ur);

atp_model_no = atp_model;
atp_model_no.lb(ur) = 0;

biomass_model_no = biomass_model;
biomass_model_no.lb(ur) = 0;

FBA = optimizeCbModel(atp_model_no) % produce 1000 atp with no nutrient because of loop

assert(FBA.f ==1000) 
FBA = optimizeCbModel(biomass_model_no);
assert (FBA.f < 1e-6) % no growth without nutrient

% remove loop 
[num,txt,raw] = xlsread('2019_01st_05_to_make_reversible.xlsx');

to_make_reversible = raw(:,1); 
assert(length(to_make_reversible)==6) % 6 reactions

[num,txt,raw] = xlsread('2019_01st_05_to_make_irreversible.xlsx');

to_make_irreversible = raw(:,1); 
assert(length(to_make_irreversible)==115) % 115 reactions

[~,iB] = ismember(to_make_irreversible,model.rxns);
[~,iA] = ismember(to_make_reversible,model.rxns);

atp_model_no.lb(iB) = 0;
atp_model_no.lb(iA) = -1000;
atp = optimizeCbModel(atp_model_no); % no atp production
assert(atp.f<1e-6)

biomass_model.lb(iB) = 0;
biomass_model.lb(iA) = -1000;
biomass_model = changeRxnBounds(biomass_model,'r_160_exchange',-100,'l');
growth = optimizeCbModel(biomass_model);
assert(growth.f>0.01) % growth 0.1022 with nitrogen -10 and glucose -20 

dataset(atp.f,growth.f,'VarNames',{'atp','growth'}) % growth 0.10221 

model = biomass_model; 