function [model_new,CNgg,lipid] = BiomassInDifConditions(model,mol_c_source,n_csource, mol_nitrogen_source,n_nitrogen_source)

%% Formulate biomass to simulate at nitrogen limited condition 
%% Find parameters 

% biomass = p protein + c carbohydrate + l lipid + 5% dna and rna 

% p + c + l = 100 -5 = 95 % 

% problem find p,c,l in different C/N ratio 


%% Input
% - model
% - mol_c_source : uptake rate of c source (for example, -10 mmol/g/h glucose)
% -n_c_source: number of carbon in the c_source. for example, glucose has 6 carbon

% Nhung Pham 2020 01st july - the script has been tested on MATLAB R2016b


%% 1. Find l. get relationship between lipid (l) in the biomass and C/N ratio from table 2 Awad et al 2019 

CNgg = [0, 12,24,30,36,48,60,90,120,180,240];
lipid = [0, 13.71, 13.48, 19.15, 15.98, 24.09, 33.78, 38.09 44.36, 41.66, 26.64];
error_bar = [0, 1.10, 1.90, 1.53, 3.4, 0.8, 4.68, 8.28, 2.56, 5.74, 3.21];

% figure (1)
% title('In-vitro lipid content at different C/N ratio')
% errorbar(CNgg, lipid,error_bar, 'bp','LineWidth',1.5)
% set(gca,'FontSize',35)
% xlabel('C/N (g/g) ')
% ylabel('% Lipid content in DCW (w/w)')

[fit, gof] = createFit(CNgg, lipid);
ylim([0 50])
xlim([0 250])
%% get the equation (1) -> get value for l 

%% Validate the equation

% data from Ykema et al 1988 Jan, Fig 2. 
CN_test = [4.636, 5.563, 7.417, 12.185, 16.689, 24.901, 46.623, 58.808, 70.993, 76.689];

lipid_invitro = [20.965, 22.355, 23.436, 19.884, 20.193, 24.517, 35.637, 40.27, 44.44, 46.911]; 

y = fit.p1 * (CN_test.^2) + fit.p2 * CN_test + fit.p3;
% 
% figure(6)
% title('In-vitro Ykema et al 1988 data validate model')
% plot (CN_test,y,'LineWidth',5)
% hold on
% sz = 100;
% scatter(CN_test,lipid_invitro,sz,'MarkerEdgeColor',[0.4660 0.6740 0.1880],...
%               'MarkerFaceColor',[0.4660 0.6740 0.1880],...
%               'LineWidth',1.5)
% ylabel('lipid content (% DCW g/g)')
% xlabel('C/N (g/g)')
% set(gca,'FontSize',30)

%% 2. Find c and p. get relationship between c and p and carbon and nitrogen (at different
% C/N )
% more carbonhydrate in the medium - more intermediate carbonhydrate in the
% cell - then it is converted to lipid -> more c in medium -> more lipid,
% intra carbonhydrate stay the same

% less nitrogen in the medium -> less protein as it is degraded to produce
% nitrogen, also not enough nitrogen to synthesize protein -> so protein is
% reduced 

% According to data in Ykema et al 1988 at CN 28 there are 11 C + 20 P + 60
% L . 60% lipid is the maximum -> the minimum C for growth is 11%. if not
% enough nitrogen for growth, extra carbon will go to lipid 
% therefore we have
c = 11; 

% we have lipid from (1) above. from there we can calculate protein. Done!

%% Formulate biomass

%% CN ratio -> lipid -> protein -> carbohydrate -> biomass reaction -> simulate 

%% Input data 


c_mol = mol_c_source / n_csource; % 1 mol glucose  = 6 mol carbon , 10 mol glucose = 10 *6 mol carbon 
c_g = c_mol * 12; % for example 10 mol glc, m carbon = n * 12 = 10 * 12 = 120 g. 

n_g = mol_nitrogen_source/n_nitrogen_source * 14; 

CNgg_input = c_g/n_g; % this is x in the (1) equation 

% find y  (lipid)
y = fit.p1 * CNgg_input * CNgg_input + fit.p2 * CNgg_input + fit.p3;

lipid = y;
% find protein 

protein = 95 - lipid - 11 ; % (5 % count for DNA, RNA and others in biomass)

% convert protein, lipid, carbohydrate percentage to coefficient for 1 g
% biomass

MW = [126, 475.14, 1465.8]; % molecular weights of protein, lipid , carbohydrate  

p_coef = (protein /MW(1)/100) * 1000; 
c_coef = (11/ MW(3)/100)*1000;
l_coef = (lipid /MW(2)/100)*1000;

%% add to model 
others = {'sulphate [cytoplasm]','riboflavin [cytoplasm]'};
[~,other_in_model] = ismember(others,model.metNames);

model_new = addReaction(model, 'CNratio_specific_Biomass', ['protein_c';'carbohydrate_c';'lipid_c'; 'dna_c';'rna_c'; model.mets(other_in_model); 'atp_c';...
    'h2o_c'; 'biomass_c';'adp_c';'h_c';'pi_c' ], [- p_coef; - c_coef; - l_coef; - 0.00918; - 0.0788; -0.000001; -0.000001; -59.276;...
    -59.276; 1; 59.276; 59.276; 59.276]);



CNgg = CNgg_input;









