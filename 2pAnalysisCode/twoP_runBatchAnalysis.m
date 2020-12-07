animal={'Plex51'};
sessionNames = {'200326','200327','200327a','200328','200331a','200401','200401a','200401b','200402'};

parfor i = 1:length(sessionNames)
    try
        twoP_scripts(animal{:},sessionNames{i});
    end
end