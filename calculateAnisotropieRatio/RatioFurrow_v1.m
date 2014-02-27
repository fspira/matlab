
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

addpath('/Users/spira/Desktop/programme/calculateAnisotropieRatio')
addpath('/Users/spira/Desktop/programme/centerOfMass')
addpath('/Users/spira/Desktop/programme/curvature')
addpath('/Users/spira/Desktop/programme/determineAngle')
addpath('/Users/spira/Desktop/programme/staging')
addpath('/Users/spira/Desktop/programme/tools')


curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell13_notch_bottomn.lsm';
tifFilenameMid = 'cell13_conv.lsm';

saveFileName = 'pos13_notch_Pole';

%load('ratioAnisoParameters.mat')

zSectionToAnalyze =6% notch casette image
zSectionMidStack = 11; % conventional image

%%%%% Load midSection

%imgMid = imread(char(tifFilenameMid));
imgMidtmp = tiffread30(char(tifFilenameMid));

imgMidtmpTmp = cat(3,imgMidtmp.data);
%imgMid = imgMidtmpTmp;

%%%%%% section required for multi stack images

imgMidtmpTmp = imgMidtmpTmp(:,:,zSectionMidStack);

voxelX = getfield(imgMidtmp,'lsm','VoxelSizeX');
voxelX_mumMid = voxelX*1000000;

for lauf = 1:3
    
    imgMid(:,:,lauf) = imgMidtmpTmp{lauf};
    
end

clear imgMidtmp
clear imgMidtmpTmp



imgOrig = tiffread30(char(tifFilename))
voxelX = getfield(imgOrig,'lsm','VoxelSizeX');
voxelX_mum = voxelX*1000000;


 truncName = findstr(tifFilename,'.lsm');
emptyFlag = isempty(truncName);

if emptyFlag ==1

 truncName = findstr(tifFilename,'.tif');

end


folderName = tifFilename(1:truncName-1);
mkdir([curdir '/' folderName]);

img = cat(3,imgOrig.data);

img = img(:,:,zSectionToAnalyze);


%imgOrigGreen = img(:,:,1);
%imgOrigRed = img(:,:,2);

[m n p] = size(img);

%p=1

%%%%% multi file tiff

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
zSectionToAnalyze = 1;
[redImg greenImg] = bkgCorrectionRedGreenExt(redImg, greenImg,zSectionToAnalyze);

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
D_Anisotropie = D;

[verticalComponent, horizontalComponent, dSum] = doVisualizeSPComponents(D,tifFilename);

%%%%%% normalize to mean itensity of each channel
%[m n p] = size(img);

S1Norm = S1;
S2Norm = S2;

for lauf = 1:p
   % S1Norm(:,:,lauf) = S1(:,:,lauf)./mean(mean(S1(:,:,lauf)));
   % S2Norm(:,:,lauf) = S2(:,:,lauf)./mean(mean(S2(:,:,lauf)));
    ratio(:,:,lauf) = S2Norm(:,:,lauf)./S1Norm(:,:,lauf);
  
end

ratio(~isfinite(ratio))=0;

%%%%%% Calculate Ratio

%Miji;
%MIJ.createImage(greenImg);

analysisFrame = zSectionToAnalyze;

imshow(imgMerge(:,:,:,analysisFrame),[])
%colormap(jet)
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
xTest =    [x-longAxis x-longAxis x+longAxis x+longAxis x-longAxis]'

yTestVer = [y-longAxis y-longAxis y+longAxis y+longAxis y-longAxis]'
xTestVer = [x-shortAxis x+shortAxis x+shortAxis x-shortAxis x-shortAxis]'


%%%%%% Total number of pixels within the box
BW_Box = (yTest(2)-yTest(1)) * (xTest(3) - xTest(2));

%%%%%% Draw the horizontal box
BWHor = roipoly(ratio(:,:,analysisFrame),xTest, yTest);

%%%%% Draw the vertical box
BWVer = roipoly(ratio(:,:,analysisFrame),xTestVer, yTestVer);


%%%%% Daw the vertical box
imgMergeVertical = doDrawROI(imgMerge,BWVer,analysisFrame);
%%%%% Daw the horizontal box
imgMergeHorizontal = doDrawROI(imgMerge,BWHor,analysisFrame);

imwrite(imgMergeVertical,'imgMergeVertical.tif','tif')
imwrite(imgMergeHorizontal,'imgMergeHorizontal.tif','tif')

%%%%% The the Roi printed into the images


BWHor = double(BWHor);
BWVer = double(BWVer);

for lauf = 1:p
    %%%%%%% Calculate D and ratio for horizontal box
    roiMaskD(:,:,lauf) = D_Anisotropie(:,:,lauf).*BWHor;
    averageRatioD(lauf) = (sum(sum(roiMaskD(:,:,lauf)))/BW_Box);
    
    roiMask(:,:,lauf) = ratio(:,:,lauf).*BWHor;
    averageRatio(lauf) = (sum(sum(roiMask(:,:,lauf)))/BW_Box);
    
    
    
    %%%% This part averages intensities within the box and subsequently
    %%%% calculates the ratio
    
   roiMask(:,:,lauf) = ratio(:,:,lauf).*BWHor;
    
   S1Mask(:,:,lauf) = S1(:,:,lauf) .*BWHor;
   S2Mask(:,:,lauf) = S2(:,:,lauf) .*BWHor;
   S1Avg(lauf) = (sum(sum(S1Mask(:,:,lauf)))/BW_Box);
   S2Avg(lauf) = (sum(sum(S2Mask(:,:,lauf)))/BW_Box);
   
   AvgRatio(lauf) = S2Avg(lauf)/S1Avg(lauf);
     
   %%%%%%% Calculate D and ratio for horizontal box
   
    roiMaskDVer(:,:,lauf) = D_Anisotropie(:,:,lauf).*BWVer;
    averageRatioDVer(lauf) = (sum(sum(roiMaskDVer(:,:,lauf)))/BW_Box);
    roiMaskVer(:,:,lauf) = ratio(:,:,lauf).*BWVer;
    averageRatioVer(lauf) = (sum(sum(roiMaskVer(:,:,lauf)))/BW_Box);
    
    
    
    %%%% This part averages first intensities within the box and
    %%%% calculates the ratio after averaging
    
   roiMaskVer(:,:,lauf) = ratio(:,:,lauf).*BWHor;
    
   S1MaskVer(:,:,lauf) = S1(:,:,lauf) .*BWVer;
   S2MaskVer(:,:,lauf) = S2(:,:,lauf) .*BWVer;
   S1AvgVer(lauf) = (sum(sum(S1MaskVer(:,:,lauf)))/BW_Box);
   S2AvgVer(lauf) = (sum(sum(S2MaskVer(:,:,lauf)))/BW_Box);
   
   AvgRatioVer(lauf) = S2AvgVer(lauf)/S1AvgVer(lauf);

end

imshow(ratio(:,:,analysisFrame))
colormap(jet)



%imwrite(dSum,[tifFilename,'_Aniso_FalseColor.tif'])
%tiffwrite_mat(dSum, [tifFilename,'_Ratio.tif'])
%tiffwrite_RBG(imgMerge,[tifFilename,'RGBImage.tif'])
tiffwrite_mat(imgMerge,[tifFilename,'RGBImage.tif']);





%%%%% High background in the mid image - therefore substract background
zSectionMidStack = 1;

%[imgMidTmp1 imgMidTmp2] = bkgCorrectionRedGreenExt(imgMid(:,:,1), imgMid(:,:,2),zSectionMidStack);

%clear imgMid
%imgMid(:,:,1) = imgMidTmp1;
%imgMid(:,:,2) = imgMidTmp2;

[distChrom distContractileRing orientation] = doAngleDistances_v1(imgMid,voxelX_mumMid,tifFilename);

stagingParameters = double([distChrom distContractileRing]);

[density flagOut predicted_time] = doMapTime(stagingParameters);

flagOut =  double(flagOut)


cd(curdir)
              saveVariables = {};

            saveVariables{1} =  distChrom;
            saveVariables{2} =  distContractileRing;
            
            
            saveVariables{3} =  orientation;
            saveVariables{4} = analysisFrame;
            saveVariables{5} = averageRatioD(analysisFrame);
            saveVariables{6} = averageRatio(analysisFrame);
            saveVariables{7} = averageRatioDVer(analysisFrame);
            saveVariables{8} = averageRatioVer(analysisFrame);
            saveVariables{9} = predicted_time;
            saveVariables{10} = flagOut;
            
            
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['chromatin distance, contractile ring distance, orientation, analysis frame, ',...
                'averageRatioD_Horizontal, averageRatio_Horizontal, averageRatioD_Vertical, averageRatio_Vertical, ',...
                'predicted_time, probability'];
            
            outid = fopen([tifFilename,'RatioAnalysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'RatioAnalysis.csv'],csvData,'roffset',1,'-append')


              saveVariables = {};

           
            saveVariables{1} = AvgRatio'; 
            saveVariables{2} = AvgRatioVer'; 
            
            
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= [' averageRatioStackHorizontal, averageRatioVertical'];
                
            
            outid = fopen([tifFilename,'RatioAnalysisStack.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'RatioAnalysisStack.csv'],csvData,'roffset',1,'-append')


  
 %   MIJ.run('closeAllWindows');
 predicted_time
 
 
   structName{1} = 'chromatin_distance'
   structName{2} =  'contractile_ring_diameter'
   structName{3} =   'orientation'
   structName{4} =  'analysis_frame'
   structName{5} = 'averageRatioD_Horizontal'
   structName{6} = 'averageRatio_Horizontal'
   structName{7} = 'averageRatioD_Vertical'
   structName{8} = 'averageRatio_Vertical' 
   structName{9} = 'predicted_time'
   structName{10} = 'probability'
   
    saveVariables = {};

            saveVariables{1} =  distChrom;
            saveVariables{2} =  distContractileRing;
            
            
            saveVariables{3} =  orientation;
            saveVariables{4} = analysisFrame;
            saveVariables{5} = averageRatioD(analysisFrame);
            saveVariables{6} = averageRatio(analysisFrame);
            saveVariables{7} = averageRatioDVer(analysisFrame);
            saveVariables{8} = averageRatioVer(analysisFrame);
            saveVariables{9} = predicted_time;
            saveVariables{10} = flagOut;
            
   
 
   
    for lauf = 1:10
        
       
        s.(saveFileName).(structName{lauf}) = struct(structName{lauf},saveVariables{lauf})
     
        
    end
    
    save([saveFileName, '.mat'],'s')%,'-append');
    %load('ratioAnisoParameters.mat')

    filename = [tifFilename, '.mat'];
    %save(filename)
