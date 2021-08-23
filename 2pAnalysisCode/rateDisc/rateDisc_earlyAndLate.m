function [cDetect, cError] = rateDisc_earlyAndLate(rInd, dInd, oInd, earlyInd, lateInd)
% short code to compute performance with and without optogenetics, before
% and after a certain point in training (turning point defined by lateDate.

% early performance - all control trials with fibers at current location
[cError(1,1,1), cError(1,1,2), cDetect(1,1)] = ...
    Behavior_wilsonError(sum(rInd & ~oInd & earlyInd), sum(dInd & ~oInd & earlyInd)); %performance with error

% early performance - all optogenetic trials with fibers at current location
[cError(1,2,1), cError(1,2,2), cDetect(1,2)] = ...
    Behavior_wilsonError(sum(rInd & oInd & earlyInd), sum(dInd & oInd & earlyInd)); %performance with error

% late performance - all control trials with fibers at current location
[cError(2,1,1), cError(2,1,2), cDetect(2,1)] = ...
    Behavior_wilsonError(sum(rInd & ~oInd & lateInd), sum(dInd & ~oInd & lateInd)); %performance with error

% late performance - all optogenetic trials with fibers at current location
[cError(2,2,1), cError(2,2,2), cDetect(2,2)] = ...
    Behavior_wilsonError(sum(rInd & oInd & lateInd), sum(dInd & oInd & lateInd)); %performance with error