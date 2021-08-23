

%% Run logistic regression on all data

twoP_lrProcessAllSessions;

%% Find empty folders

C=twoP_lrCheckResults;

%%
clear LR;

LR = twoP_lrCombineAllSessions(0);