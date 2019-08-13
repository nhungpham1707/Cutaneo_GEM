function new_model = FormulateBiomass(model)

% Nhung 26th Feb 2019
% Formulate Biomass reaction
index_biomass = strmatch('r_021_xxx',model.rxns);
biomass_component = model.metNames(find(model.S(:,index_biomass)<0));
biomass_coeff = model.S(find(model.S(:,index_biomass)<0),index_biomass);
biomass_co = full(biomass_coeff);
biomass_product = model.metNames(find(model.S(:,index_biomass)>0));

% add DNA and RNA reaction
% DNA (data for cryptococcus)

% 3.372 ATP + 0.297 dCMP + 0.297 dGMP + 0.203 dAMP + 0.203 dTMP + 3.372 H2O -> 3.372 H + 3.372 ADP + 3.372 Phosphate + DNA 

model_new = addReaction(model,'DNA_synthesis', '3.372 atp_c + 0.297 dcmp_c + 0.297 dgmp_c + 0.203 damp_c + 0.203 dtmp_c + 3.372 h2o_c -> 3.372 h_c + 3.372 adp_c + 3.372 pi_c + dna_c'); 
model_new.metNames(end) = cellstr('DNA [cytoplasm]');

% RNA (data from ecoli)
% 2.4 ATP + 0.2 CMP + 0.262 AMP + 0.216 UMP + 0.322 GMP + 2.4 H2O -> 2.4 H + 2.4 ADP + 2.4 Phosphate + RNA 
model_new = addReaction(model_new,'RNA_synthesis','2.4 atp_c + 0.2 cmp_c + 0.262 amp_c + 0.216 ump_c + 0.322 gmp_c + 2.4 h2o_c -> 2.4 h_c + 2.4 adp_c + 2.4 pi_c + rna_c'); 
model_new.metNames(end) = cellstr('RNA [cytoplasm]'); 

% Protein 

protein = {'Ala-tRNA(Ala) [cytoplasm]','Arg-tRNA(Arg) [cytoplasm]','Asn-tRNA(Asn) [cytoplasm]','Asp-tRNA(Asp) [cytoplasm]','Cys-tRNA(Cys) [cytoplasm]',...
    'Gln-tRNA(Gln) [cytoplasm]','Glu-tRNA(Glu) [cytoplasm]','Gly-tRNA(Gly) [cytoplasm]','L-glycine [cytoplasm]','His-tRNA(His) [cytoplasm]',...
    'Leu-tRNA(Leu) [cytoplasm]','Lys-tRNA(Lys) [cytoplasm]','Met-tRNA(Met) [cytoplasm]','Phe-tRNA(Phe) [cytoplasm]','Pro-tRNA(Pro) [cytoplasm]','Ser-tRNA(Ser) [cytoplasm]',...
    'Thr-tRNA(Thr) [cytoplasm]','Tyr-tRNA(Tyr) [cytoplasm]','Val-tRNA(Val) [cytoplasm]','Ile-tRNA(Ile) [cytoplasm]','Trp-tRNA(Trp) [cytoplasm]'};
protein_product = biomass_product(4:end); % stoichiometric coefficient of these metabolites also need to adjust accordingly to that of protein
[~,p_product] = ismember(protein_product,model.metNames);
p_product_coef = full(model.S(p_product,index_biomass));
p_product_coef_1_mol = p_product_coef/3.601101;

[~,p_in_model] = ismember(protein,model.metNames) ;
p_coef = full(model.S(p_in_model,index_biomass));
p_coef_for_1_mol = p_coef/3.601101;

model_new = addReaction(model_new,'Protein_synthesis',[model_new.mets(p_in_model); 'atp_c'; 'h2o_c'; model_new.mets(p_product);'protein_c'; 'adp_c'; 'pi_c'; 'h_c'],...
    [p_coef_for_1_mol; -4.306; -4.306; p_product_coef_1_mol; 1; 4.306; 4.306; 4.306]); 

% Carbonhydrate
carbon = {'(1->3)-beta-D-glucan [cytoplasm]','(1->6)-beta-D-glucan [cytoplasm]','glycogen [cytoplasm]','mannan [cytoplasm]'};
[~,c_in_model] = ismember(carbon,model.metNames);
c_coef = full(model.S(c_in_model,index_biomass));
c_coef_for_1_mol = c_coef/2.475;

model_new = addReaction(model_new,'Carbohydrate_synthesis',[model_new.mets(c_in_model);'carbohydrate_c'], [c_coef_for_1_mol;1]);

% others
others = {'sulphate [cytoplasm]','riboflavin [cytoplasm]'};
[~,other_in_model] = ismember(others,model.metNames);

%% Formulate biomass 

% model_new = addReaction(model_new, 'Biomass_nitrogen_abundant', ['protein_c';'carbohydrate_c';'lipid_c'; 'dna_c';'rna_c'; model_new.mets(other_in_model); 'atp_c';...
%     'h2o_c'; 'biomass_c';'adp_c';'h_c';'pi_c' ], [- 0.516; - 0.143; - 0.286; - 0.00567; - 0.0487; -0.000001; -0.000001; -59.276;...
%     -59.276; 1; 59.276; 59.276; 59.276]);

model_new = addReaction(model_new, 'Biomass_nitrogen_abundant', ['protein_c';'carbohydrate_c';'lipid_c'; 'dna_c';'rna_c'; model_new.mets(other_in_model); 'atp_c';...
    'h2o_c'; 'biomass_c';'adp_c';'h_c';'pi_c' ], [- 3.016; - 0.232; - 0.463; - 0.00918; - 0.0788; -0.000001; -0.000001; -59.276;...
    -59.276; 1; 59.276; 59.276; 59.276]);

% model_new = addReaction(model_new, 'Biomass_nitrogen_deletion', ['protein_c';'carbohydrate_c';'lipid_c'; 'dna_c';'rna_c'; model_new.mets(other_in_model); 'atp_c';...
%     'h2o_c'; 'biomass_c';'adp_c';'h_c';'pi_c' ], [- 0.2356; - 0.0402; - 0.677; ; - 0.0049; - 0.0422; -0.000001; -0.000001; -59.276;...
%     -59.276; 1; 59.276; 59.276; 59.276]);
model_new = addReaction(model_new, 'Biomass_nitrogen_deletion', ['protein_c';'carbohydrate_c';'lipid_c'; 'dna_c';'rna_c'; model_new.mets(other_in_model); 'atp_c';...
    'h2o_c'; 'biomass_c';'adp_c';'h_c';'pi_c' ], [- 1.587; - 0.075; - 1.263; ; - 0.00918; - 0.0788; -0.000001; -0.000001; -59.276;...
    -59.276; 1; 59.276; 59.276; 59.276]);

% Remove biomass reaction from iNL895

modelOut = removeRxns(model_new,model_new.rxns(index_biomass));

new_model = modelOut; 



