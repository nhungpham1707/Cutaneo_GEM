%% Analyze biomass reaction in the model 
%% Formulate biomass reaction in N+ and N- condition based on experimental data
%% Nhung 7th Feb 2019 


% identify biomass component in the model 
model = BuildCryptoGEM;
index_biomass = strmatch('r_021_xxx',model.rxns);
biomass_component = model.metNames(find(model.S(:,index_biomass)<0));
biomass_coeff = model.S(find(model.S(:,index_biomass)<0),index_biomass);
biomass_co = full(biomass_coeff);

% fid = fopen('2019_01_28_model_crypto_biomass_component.txt','w+');
% for i = 1: length(biomass_component)
%     fprintf(fid,'%s\t%6.4f\n',biomass_component{i},biomass_co(i))
% end

% group according to Protein / carbohydrate/ lipid/ others

protein = {'Ala-tRNA(Ala) [cytoplasm]','Arg-tRNA(Arg) [cytoplasm]','Asn-tRNA(Asn) [cytoplasm]','Asp-tRNA(Asp) [cytoplasm]','Cys-tRNA(Cys) [cytoplasm]',...
    'Gln-tRNA(Gln) [cytoplasm]','Glu-tRNA(Glu) [cytoplasm]','Gly-tRNA(Gly) [cytoplasm]','L-glycine [cytoplasm]','His-tRNA(His) [cytoplasm]',...
    'Leu-tRNA(Leu) [cytoplasm]','Lys-tRNA(Lys) [cytoplasm]','Met-tRNA(Met) [cytoplasm]','Phe-tRNA(Phe) [cytoplasm]','Pro-tRNA(Pro) [cytoplasm]','Ser-tRNA(Ser) [cytoplasm]',...
    'Thr-tRNA(Thr) [cytoplasm]','Tyr-tRNA(Tyr) [cytoplasm]','Val-tRNA(Val) [cytoplasm]','Ile-tRNA(Ile) [cytoplasm]','Trp-tRNA(Trp) [cytoplasm]'};
carbon = {'(1->3)-beta-D-glucan [cytoplasm]','(1->6)-beta-D-glucan [cytoplasm]','glycogen [cytoplasm]','mannan [cytoplasm]'};
nucleotide = {'AMP [cytoplasm]','CMP [cytoplasm]','dAMP [cytoplasm]','dCMP [cytoplasm]','dGMP [cytoplasm]','dTMP [cytoplasm]','GMP [cytoplasm]','UMP [cytoplasm]'}; 
% 'ATP [cytoplasm]' not include in the nucleotide list because its coef is
% too high compare to the rest, also it is not part of the biomass but
% rather represent the energy that cell use to grow
others = {'sulphate [cytoplasm]','riboflavin [cytoplasm]'};

total_protein = [protein,nucleotide]; 

[p_index,p_in_model] = ismember(protein,model.metNames) ;
[c,c_in_model] = ismember(carbon,model.metNames);
[n,n_in_model] = ismember(nucleotide,model.metNames);
[o,other_in_model] = ismember(others,model.metNames);

p_coef = abs(full(model.S(p_in_model,index_biomass)));
c_coef = abs(full(model.S(c_in_model,index_biomass)));
n_coef = abs(full(model.S(n_in_model,index_biomass)));
o_coef = abs(full(model.S(other_in_model,index_biomass)));

% plot 
figure (1)

pie([sum(p_coef),sum(c_coef),sum(n_coef),sum(o_coef), 1],{'Protein','Carbon','Nucleotide','Others','Lipid'})
legend('Protein','Carbon','Nucleotide','Others','Lipid')

% remove ATP to normalize the plot which count for 59.2760
n_coef2 = [n_coef(1);n_coef(3:end)]

figure(2)
p = pie([sum(p_coef),sum(c_coef),sum(n_coef2),sum(o_coef), 1],'FontSize',18)
% legend({'Protein','Carbon','Nucleotide','Others','Lipid'},'Fontsize',30)
p.FontSize = 14
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Protein','Carbon','Nucleotide','Others','Lipid'}; 
pText(1).String = strcat(txt(1),': ', percentValues(1));
pText(2).String = strcat(txt(2),': ', percentValues(2));
pText(3).String = strcat(txt(3),': ', percentValues(3));
pText(4).String = strcat(txt(4),': ', percentValues(4));
pText(5).String = strcat(txt(5),': ', percentValues(5));
set(findobj(p,'type','text'),'fontsize',30)
set(p(2:2:end),'FontSize',30)

colormap(parula)


% write to file
% fid = fopen('2019_01_29_protein_composition.txt','w+');
% for i = 1: length(protein)
%     fprintf(fid,'%s\t%6.4f\n',protein{i},p_coef(i))
% end
% 
% fid = fopen('2019_01_29_nucleotide_composition.txt','w+');
% for i = 1: length(nucleotide)
%     fprintf(fid,'%s\t%6.4f\n',nucleotide{i},n_coef(i))
% end
% 
% fid = fopen('2019_01_29_carbon_composition.txt','w+');
% for i = 1: length(carbon)
%     fprintf(fid,'%s\t%6.4f\n',carbon{i},c_coef(i))
% end
% 
% fid = fopen('2019_01_29_others_composition.txt','w+');
% for i = 1: length(others)
%     fprintf(fid,'%s\t%6.4f\n',others{i},o_coef(i))
% end


    

%% Check iMM904 model and iJO1366 model, according to Chan S et al 2017, biomass in iJO1366 MW is 1 g / mmol

index_yeast = find(iMM904.S(:,1577)<0);
yeast_biomass (:,1) = iMM904.metNames(index_yeast);
yeast_biomass (:,2) = num2cell(abs(full(iMM904.S(index_yeast,(iMM904.c~=0)))))
sum(abs(full(iMM904.S(index_yeast,(iMM904.c~=0))))) % = 124.9542 

index_ecoli = find(iJO1366.S(:,(iJO1366.c~=0))<0);
ecoli_biomass(:,1)=  iJO1366.metNames(index_ecoli);
ecoli_biomass(:,2) = num2cell(abs(full(iJO1366.S(index_ecoli,(iJO1366.c~=0)))));
sum(abs(full(iJO1366.S(index_ecoli,(iJO1366.c~=0))))) % 109.1003 

%% try to compute MW according to scripts from Chan S et al k
