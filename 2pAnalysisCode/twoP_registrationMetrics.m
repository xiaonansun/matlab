%%
% animal = 'Plex50'; session = '200322';
T = twop;
T.animal = animal; T.session = session;
T = T.s2p_output_dir;
addpath('C:\Users\Xiaonan Richard Sun\OneDrive - alumni.princeton.edu\Documents\Fiji.app\scripts'); % Adds Fiji script directory for tiff
ImageJ;
file_name = 'ops.mat';
all_ops_paths = dir(fullfile(T.dir.imagingRootDir,'**',file_name));
overwrite = false;

parfor i = 1:length(all_ops_paths)
    %%
    % i = 1;
    try
        clf
        iF = regexp(all_ops_paths(i).folder,filesep);
        animal = all_ops_paths(i).folder(iF(6)+1:iF(7)-1);
        session = all_ops_paths(i).folder(iF(8)+1:iF(9)-1);
        ops = load(fullfile(all_ops_paths(i).folder,all_ops_paths(i).name),'ops');
        ops = ops.ops{:};
        plot(ops.regDX); ax = gca;
        for j = length(ax.Children)
            ax.Children(j).LineWidth = 2;
        end
        ylabel('Shift'); xlabel('Component #');
        title([animal ' ' session ' registration metrics (ops[''regDX''])']);
        hLeg = legend('Rigid','Avg non-rigid','Max non-rigid');
        hLeg.Box = 'off';
        new_ax = fig_configAxis(ax);
        offsetAxes(new_ax);
        %%
        save_path_session = fullfile(all_ops_paths(i).folder,'regDX.png');
        saveas(gcf,save_path_session);
        save_path_analysis = fullfile(T.dir.imagingRootDir,'registration',['regDX_' animal '_' session '.png']);
        exportgraphics(gca,save_path_analysis);
        %%
        dualChan = permute(ops.regPC,[3 4 2 1]); % generate a multi-dimensional matrix compatible for viewing with imageJ
        dualChan = uint16(dualChan);
        save_path_tiff = fullfile(T.dir.imagingRootDir,'registration',['regDX_' animal '_' session '.tif']);
        %         MultiDimImg = zeros(300,400,4,5,6,'uint16');
        MultiDimImg = dualChan;
        fiji_descr = ['ImageJ=1.53t' newline ...
            'images=' num2str(size(MultiDimImg,3)*...
            size(MultiDimImg,4)) newline...
            'channels=' num2str(size(MultiDimImg,3)) newline...
            'slices=' num2str(size(MultiDimImg,4)) newline...
            'hyperstack=true' newline...
            'mode=grayscale' newline...
            'loop=false' newline...
            'min=0.0' newline...
            'max=65535.0'];  % change this to 256 if you use an 8bit image
        
        t = Tiff(save_path_tiff,'w')
        tagstruct.ImageLength = size(MultiDimImg,1);
        tagstruct.ImageWidth = size(MultiDimImg,2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsvPerSample = 16;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.Compression = Tiff.Compression.LZW;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        tagstruct.ImageDescription = fiji_descr;
        %         for frame = 1:size(MultiDimImg,5)
        for slice = 1:size(MultiDimImg,4)
            for channel = 1:size(MultiDimImg,3)
                t.setTag(tagstruct)
                t.write(im2uint16(MultiDimImg(:,:,channel,slice)));
                t.writeDirectory(); % saves a new page in the tiff file
            end
        end
        %         end
        t.close()
    end
end
