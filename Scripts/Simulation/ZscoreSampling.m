 function [Zscore,P] = ZscoreSampling(condition1,condition2);
% % 
% % Input:
% % - sampling result array for condition 1 (m reactions x n samples, sampling_array.points)
% % - sampling result array for condition 2 (m reactions x n samples, sampling_array.points)
% % - condition 1 is the reference point. compare cond 2 to 1 
% % 
% % Output: 
% % 
% % - Zscore 
% % - P probability score. if > 0.9 both changes in the same direction  [2]
% % 
% % Z = (mean1 - mean2)/ square root ((std1)^2/n1 + (std2)^2/n2)) ; n1,n2 sample size [1]


% references
% 1. Bordel, Sergio, Rasmus Agren, and Jens Nielsen. 
% ..."Sampling the solution space in genome-scale metabolic networks reveals transcriptional regulation in key enzymes." 
% ...PLoS computational biology 6.7 (2010): e1000859.
% 2. Ledesma?Amaro, Rodrigo, et al. "Genome scale metabolic modeling of the riboflavin overproducer Ashbya gossypii." 
% ...Biotechnology and bioengineering 111.6 (2014): 1191-1199.

% Nhung 25th july 2019 


sample_stats1 = calcSampleStats(condition1.points); % mean, std, mode, median, skew, kurt
sample_stats2 = calcSampleStats(condition2.points);

size_cond1 = size(condition1.points);
n_rx1 = size_cond1(1);
n_sample1 = size_cond1(2); 

size_cond2 = size(condition2.points);
n_rx2 = size_cond2(1);
n_sample2= size_cond2(2);     

for i = 1: n_rx1
Zscore (i) = (sample_stats2.mean(i) - sample_stats1.mean(i))/sqrt((sample_stats1.std(i))^2/ n_sample1 + (sample_stats2.std(i))^2/ n_sample2);
end


    