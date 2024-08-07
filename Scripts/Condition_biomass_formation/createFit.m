function [fitresult, gof] = createFit(CNgg, lipid)
%CREATEFIT(CNGG,LIPID)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : CNgg
%      Y Output: lipid
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Jun-2019 15:39:54


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( CNgg, lipid );

% Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'lipid vs. CNgg', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel CNgg
% ylabel lipid
% grid on


