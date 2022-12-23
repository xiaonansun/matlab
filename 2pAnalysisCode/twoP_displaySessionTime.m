function dt = twoP_displaySessionTime(bhv)

dt{1} = datetime(bhv.TrialStartTime(1), 'convertfrom', 'datenum', 'Format', 'MM/dd/yy HH:mm:ss.SSS');
dt{2} = datetime(bhv.TrialStartTime(end), 'convertfrom', 'datenum', 'Format', 'MM/dd/yy HH:mm:ss.SSS');

% disp('First trial started at: ');
% disp(dt{1});
% disp(['Last trial (' num2str(length(bhv.TrialStartTime)) ') started at: ' ]);
% disp(dt{2});
% disp('Session duration: ');
% disp(dt{2}-dt{1});