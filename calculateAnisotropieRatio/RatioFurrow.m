
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
addpath('/Users/spira/Documents/Matlab_scripte/tiffIO')
addpath('/Users/spira/Documents/Matlab_scripte/')
addpath('/Users/spira/Desktop/Desktop/LifeactCherry_GlGPIEgfp/131204')
addpath('/Users/spira/Documents/MATLAB_scripte/ImageProcessing/Utilities')
curdir = pwd;

 
tifFilename = 'cell4_halfStack.lsm';
tifFilenameMid = 'cell4_mid.lsm';

%%%%% Load midSection

imgMid = imread(char(tifFilenameMid));
imgMidtmp = tiffread30(char(tifFilenameMid));

voxelX = getfield(imgMidtmp,'lsm','VoxelSizeX');
voxelX_mumMid = voxelX*1000000;
clear imgMidtmp




imgOrig = tiffread30(char(tifFilename))
voxelX = getfield(imgOrig,'lsm','VoxelSizeX');
voxelX_mum = voxelX*1000000;


truncName = findstr(tifFilename,'.lsm');
folderName = tifFilename(1:truncName-1);
mkdir([curdir '/' folderName]);

img = cat(3,imgOrig.data);

[m n p] = size(img);

for lauf = 1:p
    imgOrigGreen(:,:,lauf) = img{:,1,lauf};
    imgOrigRed(:,:,lauf) = img{:,2,lauf};
end


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

%%%%% Calculate anisotropie

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

ratio(~isfinite(ratio))=0;

%%%%%% Calculate Ratio

%Miji;
%MIJ.createImage(greenImg);

analysisFrame = 4

imshow(ratio(:,:,analysisFrame))
colormap(jet)
[x,y] = ginput(1)
close all
%%%% This command draws an hourglass
%yTest = [y-10 y-10 y+10 y+10]
%xTest = [x+10 x-10 x+10 x-10];

%%%%%%%% setup of the ROI - 1.5µm x 4µm

shortAxis = round((1.5 /voxelX_mum)/2);
longAxis = round((4 / voxelX_mum)/2);

%%%%%% Edges of the horizontal box
yTest = [y-shortAxis y+shortAxis y+shortAxis y-shortAxis y-shortAxis]'
xTest = [x-longAxis x-longAxis x+longAxis x+longAxis x-longAxis]'

yTestVer = [x-longAxis x-longAxis x+longAxis x+longAxis x-longAxis]'
xTestVer = [y-shortAxis y+shortAxis y+shortAxis y-shortAxis y-shortAxis]'


%%%%%% Total number of pixels within the box
BW_Box = (yTest(2)-yTest(1)) * (xTest(3) - xTest(2));

%%%%%% Draw the horizontal box
BWHor = roipoly(ratio(:,:,analysisFrame),xTest, yTest);

%%%%% Draw the vertical box
BWVer = roipoly(ratio(:,:,analysisFrame),xTestVer, yTestVer);


%%%%% Daw the vertical box
imshow(BWVer)

BWHor = double(BWHor);
BWVer = double(BWVer);

for lauf = 1:p
    %%%%%%% Calculate D and ratio for horizontal box
    roiMaskD(:,:,lauf) = D(:,:,lauf).*BWHor;
    averageRatioD(lauf) = (sum(sum(roiMaskD(:,:,lauf)))/BW_Box);
    roiMask(:,:,lauf) = ratio(:,:,lauf).*BWHor;
    averageRatio(lauf) = (sum(sum(roiMask(:,:,lauf)))/BW_Box);
    
    
    
    %%%% This part averages first intensities within the box and
    %%%% calculates the ratio after averaging
    
    roiMask(:,:,lauf) = ratio(:,:,lauf).*BWHor;
    
    S1Mask(:,:,lauf) = S1(:,:,lauf) .*BWHor;
    S2Mask(:,:,lauf) = S2(:,:,lauf) .*BWHor;
   S1Avg(lauf) = (sum(sum(S1Mask(:,:,lauf)))/BW_Box);
   S2Avg(lauf) = (sum(sum(S2Mask(:,:,lauf)))/BW_Box);
   
   AvgRatio(lauf) = S2Avg(lauf)/S1Avg(lauf);
     %%%%%%% Calculate D and ratio for horizontal box

end

imshow(ratio(:,:,analysisFrame))
colormap(jet)



%imwrite(dSum,[tifFilename,'_Aniso_FalseColor.tif'])
tiffwrite_mat(D, [tifFilename,'_Ratio.tif'])
tiffwrite_RBG(imgMerge,[tifFilename,'RGBImage.tif'])





[distChrom distContractileRing orientation] = doAngleDistances(imgMid,voxelX_mumMid,tifFilename)
                
cd(curdir)
              saveVariables = {};

            saveVariables{1} =  distChrom;
            saveVariables{2} =  distContractileRing;
            
            
            saveVariables{3} =  orientation;
            saveVariables{4} = analysisFrame;
            saveVariables{5} = averageRatioD(analysisFrame);
            saveVariables{6} = averageRatio(analysisFrame);
                  
            
            
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['Chromatin distance, Contractile ring distance, orientation, analysis frame',...
                'AverageRatioD averageRatio'];
            
            outid = fopen([tifFilename,'RatioAnalysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'RatioAnalysis.csv'],csvData,'roffset',1,'-append')


              saveVariables = {};

           
            saveVariables{1} = averageRatioD(analysisFrame);
            saveVariables{2} = averageRatio(analysisFrame);
            saveVariables{3} = AvgRatio;    
            
            
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['AverageRatioD averageRatio FirstAverageRatio'];
                
            
            outid = fopen([tifFilename,'RatioAnalysisStack.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'RatioAnalysisStack.csv'],csvData,'roffset',1,'-append')


  
    MIJ.run('closeAllWindows');


close all

    
