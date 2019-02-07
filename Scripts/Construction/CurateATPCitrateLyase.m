function [ATP_citrate_model] = CurateATPCitrateLyase(model);
%% curate ATP-citrate lyase in model 
%model = BuildCryptoGEM;
accoa = 'accoa_c';

% there are 3 reactions that produce accoa in the cytosol in iNL895

%    'acetyl-CoA synthetase'  'acetate [cytoplasm] + ATP [cytoplasm] + coenzyme A [cytoplasm]  <=> acetyl-CoA [cytoplasm] + AMP [cytoplasm] + pyrophosphate [cytoplasm] '
%     'glycine C-acetyltransferase' 'coenzyme A [cytoplasm] + L-2-amino-3-oxobutanoate [cytoplasm]  <=> acetyl-CoA [cytoplasm] + L-glycine [cytoplasm] '
%    'hydroxymethylglutaryl CoA synthase' '3-hydroxy-3-methylglutaryl-CoA [cytoplasm] + coenzyme A [cytoplasm] + H+ [cytoplasm]  <=> acetoacetyl-CoA [cytoplasm] + acetyl-CoA [cytoplasm] + H2O [cytoplasm] '

% according to Fakas et al, 2017 there is no acetyl coa synthetase in
% oleaginous yeast. Instead accoa is produced from ATP citrate lyase 
% in iNL895, there is 'YALI0F05962g' , checking with KEGG this gene is
% indeed encode acetyl-coA synthetase. 
% in Cryptococcus curvatus annotation, there is g1308.t1 that have the Acetyl-coenzyme A synthetase N-terminus
%  there are 30 genes encode enzyme 6.2.1.1 (acetyl-coA synthetase) but
%  only g1308.t1 has confidence score of 1 . Only include this one in the
%  model for now . So c.curvatus has acetyl-coA synthetase 
% c.curvatus also has one gene for 2.3.3.8 ATP citrate lyase , g3238.t1
% confidence score 1  R00352 ADP + phosphate + acetyl-CoA + oxaloacetate = ATP +
% citrate + CoA  

model = addReaction(model,'R00352','atp_c + cit_c + coa_c -> adp_c + pi_c + aacoa_c  + oaa_c' );
model.rxnNames (end) = cellstr('ATP-citrate lyase');
model = changeGeneAssociation(model,'R00352', 'g3238.t1');
model.subsystems(end) = cellstr('Citrate cycle');

ATP_citrate_model = model; 

% add this reaction, the biomass FBA.f increase a little bit