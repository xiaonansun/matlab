        %% https://www.mathworks.com/matlabcentral/answers/389765-how-can-i-save-an-image-with-four-channels-or-more-into-an-imagej-compatible-tiff-format
        save_path_tiff = fullfile(T.dir.imagingRootDir,'registration',['regDX_' animal '_' session '.tif']);
        %         saveastiff(ops.regPC, save_path_tiff);
        dualChan = permute(ops.regPC,[3 4 2 1]); % generate a multi-dimensional matrix compatible for viewing with imageJ
        IJM.show('dualChan');
        ij.IJ.run(imp, "16-bit", ""); % convert back to 16 bit per channel
        %         IJ.run(imp, "Grays");
        imp = ij.IJ.getImage(); % need to get anew
        ij.IJ.saveAsTiff(imp, save_path_tiff);
        imp.Close()
        %         imp = ij.IJ.getImage();
        %         IJ.run(imp, "Rotate... ", "angle=90 grid=1 interpolation=Bilinear stack");
        