
clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

blue = zeros(256,3);
blue(:,3) = distvec;

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
 
tifFilename = 'Image 93_2z_Anaphase_zoom4_02.lsm';
%tifFilenameMid = 'cell6_conv.tif';

saveFileName = 'Image 93_2z_Anaphase_zoom4_02_center';

%load('ratioAnisoParameters.mat')

%zSectionToAnalyze = 1% notch casette image
%zSectionMidStack = 1; % conventional image

%%%%% Load midSection

%imgMid = imread(char(tifFilenameMid));
%imgMidtmp = tiffread30(char(tifFilenameMid));

%imgMidtmpTmp = cat(3,imgMidtmp.data);
%imgMid = imgMidtmpTmp;

%%%%%% section required for multi stack images

%timeInterval = getfield(imgOrig,'lsm','TimeOffset');
%timeInterval = timeInterval(2);
%time = {};
%time{1} = 0;

%imgMidtmpTmp = imgMidtmpTmp(:,:,zSectionMidStack);

%voxelX = getfield(imgMidtmp,'lsm','VoxelSizeX');
%voxelX_mumMid = voxelX*1000000;

%for lauf = 1:3
    
%    imgMid(:,:,lauf) = imgMidtmpTmp{lauf};
    
%end

clear imgMidtmp
clear imgMidtmpTmp



imgOrig = tiffread30(char(tifFilename))
voxelX = getfield(imgOrig,'lsm','VoxelSizeX');
voxelX_mum = voxelX*1000000;

timeInterval = getfield(imgOrig,'lsm','TimeOffset');
timeInterval = timeInterval(2);


 truncName = findstr(tifFilename,'.lsm');
emptyFlag = isempty(truncName);

if emptyFlag ==1

 truncName = findstr(tifFilename,'.tif');

end


folderName = tifFilename(1:truncName-1);
mkdir([curdir '/' folderName]);

img = cat(3,imgOrig.data);

%img = img(:,:,zSectionToAnalyze);


%imgOrigGreen = img(:,:,1);
%imgOrigRed = img(:,:,2);

[m n p] = size(img);

%p=1

%%%%% multi file tiff
planeSelector = 2;
for lauf = 1:p/2
   
    imgOrigGreen(:,:,lauf) = img{:,3,planeSelector};
    imgOrigRed(:,:,lauf) = img{:,4,planeSelector};
    imgMid(:,:,lauf) = img{:,1,planeSelector};
     planeSelector = planeSelector +2
    
end


greenImg = double(imgOrigGreen);
redImg = double(imgOrigRed);

[m n p] = size(greenImg);
%%%%% create RGB image

for lauf =1 :p

    greenStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(greenImg(:,:,lauf)),green);
    redStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(redImg(:,:,lauf)),red);
     blueStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(imgMid(:,:,lauf)),blue);
     
    imgMerge(:,:,:,lauf) = greenStackRGB(:,:,:,lauf) + redStackRGB(:,:,:,lauf)+blueStackRGB(:,:,:,lauf);
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

for lauf = 1:p
    
    [cNew1RefOut,rNew1RefOut]= doGetMijiLinescan(imgMerge(:,:,:,lauf))
    cNew1Ref{lauf} = cNew1RefOut;
    cNew1Ref{lauf} = rNew1RefOut;
    
end



for lauf = 1:p
    
    cNew1RefTmp = cNew1Ref{lauf};
    
    S1Norm(:,:,lauf) = S1(:,:,lauf)./mean(cNew1RefTmp{1});
    S2Norm(:,:,lauf) = S2(:,:,lauf)./mean(cNew1RefTmp{1});
    
    ratio(:,:,lauf) = S2(:,:,lauf)./S1(:,:,lauf);
  
end
ratio = S1;
S1 = S1Norm;
S2 = S2Norm;

%ratio(~isfinite(ratio))=0;

%%%%%% Calculate Ratio

%Miji;
%MIJ.createImage(greenImg);
for lauf = 1:p
        
        analysisFrame = zSectionToAnalyze;

        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        [Furrow1,Furrow1_S1,Furrow1_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

        xFurrow1{lauf} = x 
        yFurrow1{lauf} = y
        
        
        Furrow1_BW_BoxStore{lauf} = BW_Box;
        Furrow1_BWHor_BoxStore{lauf} = BWHor;
        Furrow1_BWVer_BoxStore{lauf} = BWVer;
        

        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        [Furrow2,Furrow2_S1,Furrow2_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

        xFurrow2{lauf} = x 
        yFurrow2{lauf} = y

        Furrow2_BW_BoxStore{lauf} = BW_Box;
        Furrow2_BWHor_BoxStore{lauf} = BWHor;
        Furrow2_BWVer_BoxStore{lauf} = BWVer;
        

        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        [Pole1, Pole1_S1, Pole1_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

        xPole1{lauf} = x 
        yPole{lauf} = y

        Pole1_BW_BoxStore{lauf} = BW_Box;
        Pole1_BWHor_BoxStore{lauf} = BWHor;
        Pole1_BWVer_BoxStore{lauf} = BWVer;

        [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        [Pole2, Pole2_S1, Pole2_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box,S1(:,:,lauf), S2(:,:,lauf))

        xPole2{lauf} = x 
        yPole2{lauf} = y

        Pole2_BW_BoxStore{lauf} = BW_Box;
        Pole2_BWHor_BoxStore{lauf} = BWHor;
        Pole2_BWVer_BoxStore{lauf} = BWVer;
    
        Furrow1Store{lauf} = Furrow1
        Furrow2Store{lauf} = Furrow2
        FurrowMeanIntensity{lauf} = (Furrow1+Furrow2)/2;

        Furrow1Raw_S1_Store(lauf) = Furrow1_S1;
        Furrow1Raw_S2_Store(lauf) = Furrow1_S2;


        Furrow2Raw_S1_Store(lauf) = Furrow2_S1;
        Furrow2Raw_S2_Store(lauf) = Furrow2_S2;

        FurrowRawMeanStore(lauf) = (Furrow2_S1 + Furrow1_S1)/2

        Pole1Raw_S1_Store(lauf) = Pole1_S1
        Pole1Raw_S2_Store(lauf) = Pole1_S2

        Pole2Raw_S1_Store(lauf) = Pole2_S1
        Pole2Raw_S2_Store(lauf) = Pole2_S2

        PoleRawMeanStore(lauf) = (Pole1_S1 + Pole1_S2)/2


        
        Pole1Store{lauf} = Pole1
        Pole2Store{lauf} = Pole2
        PoleMeanIntensity{lauf}= (Pole1+Pole2)/2;

        ratioFurrow = (Furrow1+Furrow2)/2;

        ratioPole = (Pole1+Pole2)/2;

        ratioFurrowStore{lauf} = ratioFurrow
        ratioPoleStore{lauf} = ratioPole
        
        
        contractileRingDistance{lauf} = (pdist2([xFurrow1{lauf},yFurrow1{lauf}],[xFurrow2{lauf},yFurrow2{lauf}])) * voxelX_mum


        
end

time = 1:timeInterval:p*timeInterval

for lauf = 1:p

    contractileRingDistanceTmp(lauf) = contractileRingDistance{lauf};
    ratioFurrowStoreTmp(lauf) =  ratioFurrowStore{lauf}
    ratioPoleStoreTmp(lauf) = ratioPoleStore{lauf}
    FurrowMeanIntensityTmp(lauf) = FurrowMeanIntensity{lauf}
    PoleMeanIntensityTmp(lauf) = PoleMeanIntensity{lauf}

end

plot(time,contractileRingDistanceTmp)
hold on
plot(time,ratioFurrowStoreTmp)
hold on
plot(time, ratioPoleStoreTmp)

 plot(time, FurrowMeanIntensityTmp,'-r')
 hold on
  plot(time, PoleMeanIntensityTmp,'-b')

 plot(time, FurrowRawMeanStore,'-r')
 hold on
  plot(time,PoleRawMeanStore,'-b')


      saveVariables = {};

            saveVariables{1} = time';%time{1};
            saveVariables{2} = contractileRingDistanceTmp';
            
            
            saveVariables{3} = ratioFurrowStoreTmp';
            saveVariables{4} = ratioPoleStoreTmp';
            saveVariables{5} = FurrowMeanIntensityTmp'
            saveVariables{6} = PoleMeanIntensityTmp';
            
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, contractile_ring_distance, ratio_furrow, ratio_pole',...
                'furrow_mean_intensity, pole_mean_intensity'];
            
            outid = fopen([tifFilename,'RatioAnalysis_mid.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'RatioAnalysis_mid.csv'],csvData,'roffset',1,'-append')

           save([tifFilename,'.mat'])
