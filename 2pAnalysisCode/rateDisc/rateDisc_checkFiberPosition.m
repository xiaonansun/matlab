% rateDisc_checkFiberPosition
%function to check if fiber position in optogenetic recordings is correctly
%assigned.

changeLocation = true; %if true, every location can be re-assigned. Otherwise this is just to check that locations are correct.
emptyOnly = false; %only check files without existing position
cPath = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon\';
% Animals = {'mSM80' 'mSM81' 'mSM82'};
% Animals = {'Fez7'};
% Animals = {'Plex05' 'Plex06'};
Animals = {'mSM83'};

for iAnimals = 1:length(Animals)
    
    for iChecks = 1:10 %check for files repeatedly. Sometimes the server needs a moment to be indexed correctly
        Files = ls([cPath '\' Animals{iAnimals} '\SpatialDisc\Session Data\*' Animals{iAnimals} '*']); %behavioral files in correct cPath
        if ~isempty(Files)
            break;
        end
        pause(0.1);
    end
    
    cDate = datenum(Files(:,length([Animals{iAnimals} '_SpatialDisc_'])+1:length([Animals{iAnimals} '_SpatialDisc_'])+10)); %isolate dates from Filenames
    cDate = cDate + (str2num(Files(:,length([Animals{iAnimals} '_SpatialDisc_'])+19:end-4))*0.01); %add session nr to timestamp
    [cDate,ind] = sort(cDate,'descend'); %sort counts to get the order of files to days correct. Newest file should be first in the list.
    cDate = floor(cDate);
    Files = Files(ind,:); %adjust order of filenames to get it to be chronological
    fPath = [cPath Animals{iAnimals} '\SpatialDisc\Session Data\']; %folder with behavioral data
    vPath = [cPath Animals{iAnimals} '\SpatialDisc\VideoData\']; %folder with video data
    
    % load behavioral data and check if video is available and optogenetic inactivation was done. If so check if fiber location is correct.
    for iFiles = 1:size(Files,1)
        
        load([fPath Files(iFiles,:)], 'SessionData'); %load current bhv file
        if isfield(SessionData,'Rewarded')
            SessionData.Rewarded = logical(SessionData.Rewarded);
        end
        
        useData = false;
        if isfield(SessionData, 'optoDur')
            if sum(SessionData.optoDur > 0) > 0
                useData = length(SessionData.Rewarded) > 100; % if file contains at least 100 trials;
            end
            
            if useData && emptyOnly
                useData = isempty(SessionData.Notes{1}); %ony use data if no notes are present
            end
            
            if SessionData.TrialSettings(1).SaveWebcam && useData
                
                cFile = [vPath strrep(Files(iFiles,:),'.mat','') filesep strrep(Files(iFiles,:),'.mat','_Video_0100_1.mp4')]; %first video file of face came in current recording
                if ~exist(cFile,'file')
                    cFile = [vPath strrep(Files(iFiles,:),'.mat','') filesep strrep(Files(iFiles,:),'.mat','_Video_0100_1.mj2')]; %first video file of face came in current recording
                end
                
                if exist(cFile,'file')
                    v = VideoReader(cFile);
                    cPic = readFrame(v); clear v;
                    
                    figure(99);
                    imagesc(rot90(cPic(:,:,1),3)); colormap gray; axis image
                    if ~isempty(SessionData.Notes{1})
                        title(['Current selection: ' SessionData.Notes{1}(1,:)]);
                    end
                    if changeLocation
                        SessionData.Notes{1} = questdlg('Where is the fiber?', 'Chose location', 'Frontal', 'Parietal', 'S1', 'Frontal');
                        save([fPath Files(iFiles,:)], 'SessionData'); %save notes to current bhv file
                    else
                        pause; %just pause and move to next recording when a button is pressed
                    end
                end
            end
        end
    end
end
