function [stim, actual_events, actual_duration] = get_possible_stim (stim_matrix, desired_events)

% possible_stim = stim_matrix(stim_matrix(:,1) < desired_strength + 0.05 & stim_matrix(:,1) >= desired_strength - 0.05,:);
possible_stim = stim_matrix(stim_matrix(:,5) == desired_events,:);

% if there is no possible stimulus with the dre if this is the best way to do
% it, though.esired strength
% then choose the closest one. Not su
if isempty(possible_stim)
    [min_difference, ind] = min(abs(stim_matrix(:,5) - desired_events));
    this_stim = stim_matrix(ind,:);
    % if there are possible stimuli, then pick one randomly
else
    ind = randperm(size(possible_stim,1));
    this_stim = possible_stim(ind(1),:);
end

% create the corresponding stimulus vector
stim = [repmat(1, 1, this_stim(2)), repmat(2, 1, this_stim(3))];
stim
rand_indexes = randperm(length(stim));

stim = stim(rand_indexes);
actual_events = this_stim(5);
actual_duration = this_stim(4);

end