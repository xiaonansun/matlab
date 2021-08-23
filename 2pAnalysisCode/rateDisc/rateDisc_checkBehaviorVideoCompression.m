recs = dir(pwd);
recs(1).isdir = false;
recs(2).isdir = false;
a = cat(1,recs([recs.isdir]).name);
b=false(1, size(a,1));
for x = 1 : size(a,1)
cFiles = dir([pwd filesep a(x,:)]);
for y = 1 : length(cFiles)
if strcmpi(cFiles(y).name, 'motionbhvOpts.mat')
b(x) = true;
end
end
end
disp(a(~b,:));
