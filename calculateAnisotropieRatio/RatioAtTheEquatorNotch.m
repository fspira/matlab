%%%%%% CaluclateAnisotropie at the cell equator. This script will first
%%%%%% segment the cell - then allow the user to validate segmentation.
%%%%%% Second user selection of centrall furrow region as well as a region
%%%%%% with simillar curvatore outside the central furrow. In addition
%%%%%% chromatin-chromatin distance hast to be selected by the user Third, curvature
%%%%%% is computed and furrow vs. non furrow curvature is computed - to
%%%%%% assure comparable curvature. Fourth, computed furrow ingression
%%%%%% diameter and chromatin-chromatin distance.


clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

javaaddpath '/Applications/Fiji.app/scripts/Miji.m'
javaaddpath '/Applications/Fiji.app/jars/ij-1.47b.jar'
addpath('/Applications/Fiji.app/scripts')

addpath('/Applications/Fiji.app/scripts')
addpath('/Users/spira/Documents/Matlab_scripte/Image_Processing_utils')
addpath('/Users/spira/Desktop/programme')
addpath('/Users/spira/Documents/Matlab_scripte/tiffIO')
addpath('/Users/spira/Documents/Matlab_scripte/')
addpath('/Users/spira/Desktop/Desktop/LifeactCherry_GlGPIEgfp/131204')
addpath('/Users/spira/Documents/MATLAB_scripte/ImageProcessing/Utilities')
curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell12_FullStack.lsm';
tifFilenameMid = 'cell12_mid.lsm';

zSectionToAnalyze = 19% notch casette image
zSectionMidStack = 1; % conventional image

%%%%% Load midSection

imgMid = imread(char(tifFilenameMid));
imgMidtmp = tiffread30(char(tifFilenameMid));

imgMidtmpTmp = cat(3,imgMidtmp.data);
imgMid = imgMidtmpTmp;
imgMidtmpTmp = imgMidtmpTmp(:,:,zSectionMidStack);

voxelX = getfield(imgMidtmp,'lsm','VoxelSizeX');
voxelX_mumMid = voxelX*1000000;

%for lauf = 1:3
    
%    imgMid(:,:,lauf) = imgMidtmpTmp{lauf};
    
%end

clear imgMidtmp
clear imgMidtmpTmp



imgOrig = tiffread30(char(tifFilename))
%voxelX = getfield(imgOrig,'lsm','VoxelSizeX');
%voxelX_mum = voxelX*1000000;


 truncName = findstr(tifFilename,'.lsm');
emptyFlag = isempty(truncName);

if emptyFlag ==1

 truncName = findstr(tifFilename,'.tif');

end


folderName = tifFilename(1:truncName-1);
mkdir([curdir '/' folderName]);

img = cat(3,imgOrig.data);
img = img(:,:,zSectionToAnalyze);
imgOrigGreen = img{1};
imgOrigRed = img{2};

[m n p] = size(img);

p=1

%for lauf = 1:p
%    imgOrigGreen(:,:,lauf) = img{:,1,lauf};
%    imgOrigRed(:,:,lauf) = img{:,2,lauf};
%end


greenImg = double(imgOrigGreen);
redImg = double(imgOrigRed);

%%%%% create RGB image

for lauf =1 :p

    greenStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(greenImg(:,:,lauf)),green);
    redStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(redImg(:,:,lauf)),red);
             
    imgMerge(:,:,:,lauf) = greenStackRGB(:,:,:,lauf) + redStackRGB(:,:,:,lauf);
end

[redImg greenImg] = bkgCorrectionRedGreen(redImg, greenImg);

%%%%% Setting negative values to small values, to avoid division through zero
greenNormZero = greenImg > 0;
greenNorm = double((greenImg .* greenNormZero));

redNormZero = redImg > 0;
redNorm = double(redImg .* redNormZero);

%%%%% change directory
cd ([curdir '/' folderName]);
workdir = pwd;


S1 = greenNorm;
S2 = redNorm;

D=((S1-1.079.*S2)./(S1+(2*1.079).*S2)).*255;
D(~isfinite(D)) = 0;
[D, dNegInv, dSum] = doVisualizeSPComponents(D,tifFilename);

%%%%%% normalize to mean itensity of each channel
[m n p] = size(img);
for lauf = 1:p
    S1Norm(:,:,lauf) = S1(:,:,lauf)./mean(mean(S1(:,:,lauf)));
    S2Norm(:,:,lauf) = S2(:,:,lauf)./mean(mean(S2(:,:,lauf)));
    ratio(:,:,lauf) = S2Norm(:,:,lauf)./S1Norm(:,:,lauf);
  
    
end

