function [fdata, cor, type] = myTrendPol2(Data,n)
% this function does the trendelimination
%
% [fdata cor type] = fitcurve(Data,n)
%
% it reutrns the data without trend in the form:
%       first column time
%       all other collums data
% the correlation-coeffiant of the trend elimination (cor)
% and a variable type containing information abaut trendelimnation
%
% it needs data (Data) in the form:
%       first column time
%       all other columns data
% it needs the number of the datarow which is trend eliminated (n)
% 
%    written by Daniel Lohse

t = Data(:,1);
data = Data(:,n);

% Fit mit quadratischer Funktion
coef = polyfit(t,data,2);
trend = coef(1)*t.^2 + coef(2)*t + coef(3);
r = corrcoef(data,trend);
cor = r(1,2);

fdata(:,1) = t;
% die Residuien werden zurueck gegeben
fdata(:,2) = data - trend;

if cor <= 0.001
  type = -1;
  fdata(:,2) = data;
else
  type = 0;
end