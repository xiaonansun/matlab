function [timeStamps,data] = Behavior_LoadVideo(cPath)
% short routine to load data from webcams in behavioral setups.
% cPath is the path of the file that should be opened.

fID = fopen(cPath);
hSize = fread(fID,1,'double'); %header size
timeStamps = fread(fID,hSize,'double'); %Metadata. Default is 1:hSize = Absolute timestamps for each frame
dSize = fread(fID,1,'double'); %numer of data array dimensions
dSize = fread(fID,dSize,'double'); %data size
data = fread(fID,[prod(dSize),1],'uint8=>uint8'); %get data. Last 4 header values should contain the size of the data array.
data = reshape(data,dSize'); %reshape data into matrix
fclose(fID);
        
end