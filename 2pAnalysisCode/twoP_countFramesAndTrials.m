function twoP_countFramesAndTrials
%%
bin_filepath = "Y:\data\richard_s2p_npy\CSP27\imaging\20200319a\200319_005.bin";

if iscell(bin_filepath) && length(bin_filepath) >= 2
[volt, ~] = readMOMAnalog(bin_filepath);
        slowGalvo{j} = volt(slowGalvoCh, :);
        [frameStarts{j}, incompleteFrames{j}] = parseSlowGalvo(slowGalvo{j});
        numOfGalvoFrames(j) = length(frameStarts{j});
        numOfIncompleteFrames(j) = length(incompleteFrames{j});
else
    [volt, ~] = readMOMAnalog(bin_filepath);
        slowGalvo = volt(slowGalvoCh, :);
        [frameStarts, incompleteFrames] = parseSlowGalvo(slowGalvo);
        numOfGalvoFrames = length(frameStarts);
        numOfIncompleteFrames = length(incompleteFrames);
end

disp(['Number of galvo frames: ' num2str(numOfGalvoFrames)]);

%%
tif_filepath = "\\grid-hs\churchland_nlsas_data\data\CSP27\imaging\20200325\200325_001_044.TIF";
tif_info = imfinfo(tif_filepath);
tif_frame_count = length({info.Width})/2;
disp(['Number of frames in file: ' num2str(tif_frame_count)])

%% 
spks_filepath ="Y:\data\richard_s2p_npy\CSP27\imaging\20200319a\suite2p\plane0\spks.npy";
spks = readNPY(spks_filepath);
disp(['Number of frames in spks: ' num2str(length(spks))]); 