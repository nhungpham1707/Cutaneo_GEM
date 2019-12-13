function model = ConvertID(model);

%% this function convert model in-house ID to ID in standard database (BiGG, CheBI, and KEGG)
%% Input
%- model

%% output

%-model with new standard ID 

%% Nhung 17-08-2018 
% model = BuildFinalCryptoModel; 
ori_dir = pwd;
% fileName = '15_08_2018_Map_CheBI_BiGG.xlsx';
fileName = '01st_08_2019_Map_CheBI_BiGG.xlsx'; % Nhung changed some CheBI ID to BiGG , i.e. citrate, oxaloacetate

% load data
dir = fileparts(which(fileName));
cd (dir) 
[NUM,TXT,RAW]=xlsread(fileName,1); % 15_08_2018_Map_CheBI_BiGG

data = RAW(2:end,:); % remove header met name in model, curate ID, comment, met name in standard databases, IDs in standard databases

% prepare model data
w3 = regexprep(model.metNames,'\[acyl-carrier protein\]','(acyl-carrier protein)','ignorecase');
w4(:,1) = model.mets;
w4(:,2) = model.metNames;
w4(:,3) = regexp(w3, '\[[a-zA-z]+\s*[a-zA-z]*\]','match', 'once'); % take compartments


compartment (:,1) = {'[cell envelope]','[cytoplasm]','[extracellular]','[mitochondrion]','[nucleus]','[peroxisome]','[endoplasmic reticulum]',...
    '[vacuole]','[Golgi]','[lipid particle]','[vacuolar membrane]'}; 
compartment (:,2) = {'ce','c','e','m','n','x','r','v','g','l','vm'}; % ce for cell evelop, l for lipid particle, there are no standard format

for i = 1: length(compartment)
    index = strmatch(compartment(i,1),w4(:,3));
    w4(index,4) = compartment(i,2);
end

w4(:,5) = regexprep(w3,'\[[a-zA-z]+\s*[a-zA-z]*\]','','ignorecase'); % remove compartment
   
% map data to model

for i = 1: length(data)
    index = strmatch(data(i,1),w4(:,5));
    if isempty(find(strcmp(data(i,3),'c')))
        w4(index,6) = cellstr(num2str(data{i,5}));
    elseif ~ isempty(find(strcmp(data(i,3),'c')))
        w4(index,6) = cellstr(num2str(data{i,2}));
    end
end

for i = 1: length(w4);
    w4(i,7) = strcat(w4(i,6),'_',w4(i,4));
end

% change met IDs in model
ori_model = model;
for i = 1: length(model.mets)
    if isempty(strmatch(regexprep(w4(i,7),'_','','ignorecase'),w4(i,4),'exact')); % met ID that did not match standard database
        model.mets(i) =  w4(i,7);
    else
        model.mets(i) = strcat(model.mets(i),w4(i,7));
    end
end

model = model; 

cd (ori_dir)        
