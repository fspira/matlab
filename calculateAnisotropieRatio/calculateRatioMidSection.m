
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

load('voxelX_mum.mat');
load('voxelX_mumMid');
 
tifFilename = 'cell6_notch_mid.tif';
tifFilenameMid = 'cell6_conv.tif';

saveFileName = 'cell6_notch_Mid_center';

%load('ratioAnisoParameters.mat')

zSectionToAnalyze = 1% notch casette image
zSectionMidStack = 1; % conventional image

%%%%% Load midSection

%imgMid = imread(char(tifFilenameMid));
imgMidtmp = tiffread30(char(tifFilenameMid));

imgMidtmpTmp = cat(3,imgMidtmp.data);
imgMid = imgMidtmpTmp;

%%%%%% section required for multi stack images

%timeInterval = getfield(imgOrig,'lsm','TimeOffset');
%timeInterval = timeInterval(2);

time = 0;

%imgMidtmpTmp = imgMidtmpTmp(:,:,zSectionMidStack);

%voxelX = getfield(imgMidtmp,'lsm','VoxelSizeX');
%voxelX_mumMid = voxelX*1000000;

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

%img = img(:,:,zSectionToAnalyze);


imgOrigGreen = img(:,:,1);
imgOrigRed = img(:,:,2);

[m n p] = size(img);

p=1

%%%%% multi file tiff

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

%Miji

[cNew1Ref,rNew1Ref]= doGetMijiLinescan(imgMerge)



for lauf = 1:p
    S1Norm(:,:,lauf) = S1(:,:,lauf)./mean(cNew1Ref{lauf});
    S2Norm(:,:,lauf) = S2(:,:,lauf)./mean(cNew1Ref{lauf});
    
    ratio(:,:,lauf) = S2(:,:,lauf)./S1(:,:,lauf);
  
end

S1 = S1Norm;
S2 = S2Norm;

%ratio(~isfinite(ratio))=0;

%%%%%% Calculate Ratio

%Miji;
%MIJ.createImage(greenImg);
for lauf = 1:p
        
        analysisFrame = zSectionToAnalyze;

        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge,analysisFrame,ratio,voxelX_mum);
        [Furrow1] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1, S2)

        xFurrow1{lauf} = x 
        yFurrow1{lauf} = y

        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge,analysisFrame,ratio,voxelX_mum);
        [Furrow2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1, S2)

        xFurrow2{lauf} = x 
        yFurrow2{lauf} = y


        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge,analysisFrame,ratio,voxelX_mum);
        [Pole1] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1, S2)

        xPole1{lauf} = x 
        yPole{lauf} = y


        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge,analysisFrame,ratio,voxelX_mum);
        [Pole2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1, S2)

        xPole2{lauf} = x 
        yPole2{lauf} = y

    
        Furrow1Store{lauf} = Furrow1
        Furrow2Store{laauf} = Furrow2
        FurrowMeanIntensity{lauf} = (Furrow1+Furrow2)/2;
        
        Pole1Store{lauf} = Pole11
        Pole2Store{lauf} = Pole2
        PoleMeanIntensity{lauf}= (Pole1+Pole2)/2;

        ratioFurrow = (Furrow1+Furrow2)/2;

        ratioPole = (Pole1+Pole2)/2;

        ratioFurrowStore{lauf} = ratioFurrow
        ratioPoleStore{lauf} = ratioPole
        
        
        contractileRingDistance{lauf} = (pdist2([xFurrow1{lauf},yFurrow1{lauf}],[xFurrow2{lauf},yFurrow2{lauf}])) * voxelX_mum


        
end

      saveVariables = {};

            saveVariables{1} = time;
            saveVariables{2} = contractileRingDistance;
            
            
            saveVariables{3} = ratioFurrowStore;
            saveVariables{4} = ratioPoleStore;
            saveVariables{5} = FurrowMeanIntensity
            saveVariables{6} = PoleMeanIntensity
            
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, contractile ring distance, ratio furrow, ratio pole, ',...
                'furrow mean intensity, pole mean intensity'];
            
            outid = fopen([tifFilename,'RatioAnalysis_mid.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'RatioAnalysis_mid.csv'],csvData,'roffset',1,'-append')

           save(tifFilename)
