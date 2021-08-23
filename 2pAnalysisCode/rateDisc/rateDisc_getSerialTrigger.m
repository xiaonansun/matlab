
sync = single(sync);
numChans = 16;
targChan = 16;

dataOut = false(length(sync),1);
for x = 1 : length(sync)
    cData = d2b((sync(x)), 16);
    dataOut(x) = cData(targChan);
    if rem(x,length(sync) / 100) == 0
        disp(x / length(sync));
    end
end

% isolate digital lines from sync vector
[f,e]=log2(max(sync));
a = logical(rem(floor(sync*pow2(1-max(numChans,e):0)),2));
dataOut = a(:,targChan); %get target channel


% get events
trigOn = find(diff(dataOut) == -1); %all times lines switches to true
newTrig = trigOn([1; find(diff(trigOn) > std(diff(trigOn)) * 2) + 1; length(trigOn)]); %onset of new stim trigger
for x = 1 : length(newTrig) - 1
    
    plot(dataOut(newTrig(x) : newTrig(x) + 100)); ylim([0 1.2]);

    cTrigs = find(diff(dataOut(newTrig(x) : newTrig(x+1))) ~= 0); %times on line change
    cTrigs = round(diff(cTrigs) ./ (sRate / targRate)); %nr of bits in different states. First bit is true.
    
    % get first character
    a = []; checker = 1; Cnt = 0;
    while length(a) < 9
        Cnt = Cnt + 1;
        a = [a repmat(num2str(checker), 1, cTrigs(Cnt))];
        checker = double(~logical(checker));
    end
    
    temp(x) = (char(bin2dec(a(2:end))));
    temp1(x) = sum(dataOut(newTrig(x) : newTrig(x) + 1200));
    
%    pause;
    
end

checker = true;
while checker
    
    
cTrig = find(dataOut < 1, 1) - 1; %trigger


temp(:,x) = dataOut(lastStop + cTime+1 : lastStop + cTime + 1200);
lastStop = cTime + 12000;

end

%
c = (diff([1; temp(:,6)] == 0)) ~= 0;

diff(find(c)) ./ (sRate / targRate)

frames = 

