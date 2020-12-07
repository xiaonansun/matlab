function M = create_stim_matrix(shorti, longi, event, duration)

% let's go through all the possible stimulus for these settings
% by iterating the number of LONG intervals

max_nr_shorti = floor((duration - event) / (shorti + event));

nr_shorti = [0 : max_nr_shorti];
nr_longi  = floor(((duration - event) -  nr_shorti .* (shorti + event)) ./ (longi + event));

stim_strength = nr_shorti ./ (nr_shorti + nr_longi);
actual_duration = nr_longi * (longi + event) + nr_shorti * (shorti + event) + event;

M = [stim_strength', nr_shorti', nr_longi', actual_duration', nr_longi' + nr_shorti' + 1];

% take only stimuli that are longer than 930 ms
M = M(M(:,4) > 930,:);

% M = unique(M, 'rows');

end