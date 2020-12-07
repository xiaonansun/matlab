%% merging 2p stacks
clear all;
i1Path = "H:\Fez51\200331_017\200331_017_002_001.TIF";
i2Path = "H:\Fez51\200331_017\200331_017_002_002.tif";
i3Path = "H:\Fez51\200331_017\200331_017_002_003.TIF";
x_dim=512; y_dim=512;

info1 = imfinfo(i1Path);
num_images1 = numel(info1);
stack1 = zeros(x_dim,y_dim,num_images1);
parfor k = 1:num_images1
    stack1(:,:,k) = imread(i1Path, k);
end

info2 = imfinfo(i2Path);
num_images2 = numel(info2);
stack2 = zeros(x_dim,y_dim,num_images2);
parfor k = 1:num_images2
    stack2(:,:,k) = imread(i2Path, k);
end

info3 = imfinfo(i3Path);
num_images3 = numel(info3);
stack3 = zeros(x_dim,y_dim,num_images3);
parfor k = 1:num_images3
    stack3(:,:,k) = imread(i3Path, k);
end
%%
stack002=cat(3,stack1,stack2,stack3);
stack002=uint16(stack002);
stackObj= Tiff('H:\Fez51\200331_017\200331_017_002a.TIF','w');
write(stackObj,stack002);
saveastiff(stack002,'H:\Fez51\200331_017\200331_017_002b.TIF');
%%

i1 = imread(i1Path);
i2 = imread(i2Path);
i3 = imread(i3Path);

obj1 = Tiff(i1Path,'r'); i1 = read(obj1);
obj2 = Tiff(i2Path,'r'); i2 = read(obj2);
obj3 = Tiff(i3Path,'r'); i3 = read(obj3);

