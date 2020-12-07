function alldata = getStableTrialsForUnit(alldata, unit)
% alldata = getStableTrialsForUnit(alldata, unit)
%
% This function takes an alldata struct that has already had stability of
% the units assessed with markAlldataInstability(). It returns an alldata
% struct containing only the trials that were stable for unit 'unit'.

try
  stable = arrayfun(@(x) x.units(unit).stable, alldata);
catch
  warning('getStableTrialsForUnit:needStabilityAssessed', ...
    'Must run markAlldataInstability() on the data before using this function');
end

alldata = alldata(stable);
