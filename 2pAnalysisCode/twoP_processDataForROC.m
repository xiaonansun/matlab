tr = 1/30.9;
firstSideTryAl = struct;
a = permute(data.neural,[2 1 3]);
aNorm = bsxfun(@rdivide,abs(a),sqrt(sum(a.*a,2)));
firstSideTryAl.traces = aNorm;
firstSideTryAl.eventI = data.trialStimFrame;
firstSideTryAl.time = (1-data.trialStimFrame)*tr:tr:(size(data.neural,2)-data.trialStimFrame)*tr;

