function new_model = AssignSubsystem(model); 
%% Asign subsystem for reactions that may influence lipid synthesis 
%% ref. KEGG pathway

%% Nhung 14th Aug 2019 

% replace [] with empty string 

model.subSystems(cellfun(@isempty,model.subSystems)) = cellstr('');

% Find reactions in aa metabolism (based on iMM904 model)

aa = {'S_Alanine_and_Aspartate_Metabolism','S_Arginine_and_Proline_Metabolism',...
'S_Asparagine_metabolism', 'S_Cysteine_Metabolism','S_Glutamate_metabolism',... 
'S_Glutamine_Metabolism', 'S_Glycine_and_Serine_Metabolism', 'S_Histidine_Metabolism',...
'S_Methionine_Metabolism', 'S_Threonine_and_Lysine_Metabolism',...
'S_Tyrosine__Tryptophan__and_Phenylalanine_Metabolism', ...
'S_Valine__Leucine__and_Isoleucine_Metabolism'}; 

ala = {'r_0229','r_0233','r_0235','r_0236','r_0237'}; % out of 9 rx in iMM904

arg = {'r_0130','r_0731','r_0217','r_0224','r_0225','r_0226','r_0277','r_0654',...
    'r_0655','r_0656','r_0661','r_0662','r_0659','r_0660','r_0789','r_0790','r_0792','r_0791',...
    'r_0936','r_0657'}; % out of 33
asp = {'r_0229'}; % out of 2 (don't have the one in extracellular space)
cys = {'r_0024' ,'r_0974'}; % out of 10
glu = {'r_0729','r_0296','r_0297','r_0516','r_0508', 'r_0510','r_1004' }; % out of 17
gln = {'r_0221','r_0516','r_0515','r_0513'}; % out of 4
gly = {'r_0174','r_0338','r_0541','r_0546','r_0543','r_0544','r_0545','r_0539',...
    'r_0540', 'r_0538',  'r_0586','r_0587','r_0588','r_0893','r_0664'}; % out of 19
his = {'r_0208','r_0580','r_0575','r_0576','r_0577','r_0604','r_0605','r_0881','r_0882',...
    'r_0008'}; % out of 14
met = { 'r_0160','r_0159','r_0782','r_0783', 'r_0339','r_0340','r_0589','r_0702',...
    'r_0584', 'r_0786'}; % out of 20
threo ={'r_0651','r_0970', 'r_0645','r_0649','r_0537','r_0581','r_0585','r_0033','r_0668',...
    'r_0969','r_0648','r_0667','r_1026','r_1027'}; % out of 19
tyr = {'r_0187','r_0189','r_0203','r_0195','r_0204'};
val = {};


% assign subsystem to reactions in model 

model.subSystems(ismember(model.rxns,ala)) = cellstr(aa(1));
model.subSystems(ismember(model.rxns,arg)) = cellstr(aa(2));
model.subSystems(ismember(model.rxns,asp)) = cellstr(aa(3));
model.subSystems(ismember(model.rxns,cys)) = cellstr(aa(4));
model.subSystems(ismember(model.rxns,glu)) = cellstr(aa(5));
model.subSystems(ismember(model.rxns,gln)) = cellstr(aa(6));
model.subSystems(ismember(model.rxns,gly)) = cellstr(aa(7));
model.subSystems(ismember(model.rxns,his)) = cellstr(aa(8));
model.subSystems(ismember(model.rxns,met)) = cellstr(aa(9));
model.subSystems(ismember(model.rxns,threo)) = cellstr(aa(10));
model.subSystems(ismember(model.rxns,tyr)) = cellstr(aa(11));
% model.subSystems(ismember(model.rxns,val)) = cellstr(aa(12));



