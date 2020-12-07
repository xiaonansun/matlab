function [MA, MS, dupes] = deduplicateM(MA, MS, rThresh)

tetrodes = [MA.tetrode];
uTetrodes = unique(tetrodes);

files = {MA.filename};
uFiles = unique(files);

dupes = zeros(1, length(MA));

% Slide along days, compare within and adjacent days
for f = 1:length(uFiles) - 1
  for t = uTetrodes
    
    % Find units in either this day or the next one, on this tetrode
    theseUnitsToday = find(strcmp(files, uFiles{f}) & tetrodes == t);
    theseUnitsTomorrow = find(strcmp(files, uFiles{f+1}) & tetrodes == t);
    theseUnits = [theseUnitsToday theseUnitsTomorrow];
    
    % Get these units from M
    thisM = MA(theseUnits);
    
    if isempty(thisM), continue; end;
    
    % Make the A matrix from only these units
    A = makeAMatrixFromM(thisM, 0);
    
    % Get rid of NaN rows from missing conditions
    A = A(~any(isnan(A), 2), :);
    
    % Corr matrix
    R = corr(A);
    
    % Find R^2's that exceed threshold. Autocorrs ignored later.
    aboveThresh = (R >= rThresh);
    
    % Ignore within-day pairs for tomorrow (so we don't get them twice)
    nToday = length(theseUnitsToday);
    nTomorrow = length(theseUnitsTomorrow);
    aboveThresh(nToday+1:end, nToday+1:end) = false(nTomorrow);
    
    theseDupes = zeros(size(R));
    for i = 1:length(thisM)-1
      possibleDupes = i + find(aboveThresh(i, i+1:end));
      
      if ~isempty(possibleDupes)
        for p = possibleDupes
          
          hf1 = decisionPSTH(thisM(i), 1);
          addCenterText(sprintf('R^2 = %0.3f', R(i, p)));
          hf2 = decisionPSTH(thisM(p), 1);
          addCenterText(sprintf('R^2 = %0.3f', R(i, p)));
          
          isDupe = input(sprintf('%d possible pairs with this unit. Duplicates? (0 for No, 1 for Yes): ', ...
            length(possibleDupes)));
          
          while isempty(isDupe) || ~isnumeric(isDupe) || ~ismember(isDupe, [0 1])
            fprintf('Must enter a 0 or 1\n');
            isDupe = input(sprintf('%d possible pairs with this unit. Duplicates? (0 for No, 1 for Yes): ', ...
              length(possibleDupes)));
          end
          
          % Mark the lower rated unit as the dupe, or lower SNR if a tie
          if isDupe
            igtp = thisM(i).rating >= 3 && thisM(p).rating <= 3;
            pgti = thisM(i).rating <= 3 && thisM(p).rating >= 3;
            if igtp || ~pgti && thisM(i).SNR >= thisM(p).SNR
              theseDupes(i, p) = isDupe;
            else
              theseDupes(p, i) = isDupe;
            end
          end
          
          try close(hf1); end
          try close(hf2); end
        end
      end
    end
    
    dupes(theseUnits) = dupes(theseUnits) | (sum(theseDupes) > 0);
    
  end
end

% Cull the duplicates
dupes = find(dupes);
MA(dupes) = [];
MS(dupes) = [];


function addCenterText(str)

xLims = get(gca, 'XLim');
yLims = get(gca, 'YLim');

text(mean(xLims), yLims(2) * 0.9, str, 'HorizontalAlignment', 'center');

