
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


anaOnset =4;
ingressionFrame = 22;
lastFrameToConsider = 26;

curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell5_H2B_mRFP_SiRActin_Blebbistatin75muM_vertical.lsm';
%tifFilenameMid = 'cell6_conv.tif';

saveFileName = 'cell5_H2B_mRFP_SiRActin_Blebbistatin75muM_CR';

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
planeSelector = 1;


for lauf = 1:p
    
    

    imgOrigGreenTmp = img(:,:,(lauf));
    imgOrigGreen(:,:,lauf) = imgOrigGreenTmp{3};
    
    
    imgOrigRedTmp = img(:,:,(lauf));
    imgOrigRed(:,:,lauf) = imgOrigRedTmp{4};
   
     
    imgMidTmp = img(:,:,(lauf));
    imgMid(:,:,lauf) = imgMidTmp{1};
    
    %planeSelector = planeSelector +2
    
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

greenNorm = greenImg;
redNorm = redImg;

%%% Because floating crystals in the background not measured
%[redImg greenImg] = bkgCorrectionRedGreenExt(redImg, greenImg,zSectionToAnalyze);

%%%%% Setting negative values to small values, to avoid division through zero
%greenNormZero = greenImg > 0;
%greenNorm = double((greenImg .* uint32(greenNormZero)));

%redNormZero = redImg > 0;
%redNorm = double(redImg .* uint32(redNormZero));

%%%%% change directory
cd ([curdir '/' folderName]);
workdir = pwd;

%%%%% Calculate anisotropie

S1 = greenNorm;
S2 = redNorm;

%[verticalComponent, horizontalComponent, dSum] = doVisualizeSPComponents(D,tifFilename);


%%%%%% normalize to mean itensity of each channel
%[m n p] = size(img);

S1Norm = S1;
S2Norm = S2;
p = lastFrameToConsider;

%Miji
% 
% for lauf = 1:p
%     
%     [cNew1RefOut,rNew1RefOut]= doGetMijiLinescan(imgMerge(:,:,:,lauf))
%     cNew1Ref{lauf} = cNew1RefOut;
%     cNew1Ref{lauf} = rNew1RefOut;
%     
% end
% 
% 
% 
% for lauf = 1:p
%     
%     cNew1RefTmp = cNew1Ref{lauf};
%     
%     S1Norm(:,:,lauf) = S1(:,:,lauf)./mean(cNew1RefTmp{1});
%     S2Norm(:,:,lauf) = S2(:,:,lauf)./mean(cNew1RefTmp{1});
%     
%     ratio(:,:,lauf) = S2(:,:,lauf)./S1(:,:,lauf);
%   
% end
ratio = S1;
S1 = S1Norm;
S2 = S2Norm;

%ratio(~isfinite(ratio))=0;

%%%%%% Calculate Ratio

%Miji;
%MIJ.createImage(greenImg);
[m n k p] = size(imgMerge);

for lauf = 7:p
        
        analysisFrame = zSectionToAnalyze;

        [BW_Box_Out BWHor_Out BWVer_Out x_Out y_Out] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        
        BW_Box = BW_Box_Out{1};
        BWHor = BWHor_Out{1};
        BWVer  =BWVer_Out{1};
        
        x = x_Out{1};
        y = y_Out{1};
        
        xFurrow1{lauf} = x;
        yFurrow1{lauf} = y;
        
        
        Furrow1_BW_BoxStore{lauf} = BW_Box;
        Furrow1_BWHor_BoxStore{lauf} = BWHor;
        Furrow1_BWVer_BoxStore{lauf} = BWVer;
        
        
        [Furrow1,Furrow1_S1,Furrow1_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box, S1(:,:,lauf), S2(:,:,lauf))


       % [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        
        BW_Box = BW_Box_Out{2};
        BWHor = BWHor_Out{2};
        BWVer  =BWVer_Out{2};
        
        x = x_Out{2};
        y = y_Out{2};
       
       [Furrow2,Furrow2_S1,Furrow2_S2] = doCalculateRatioBoundingBox(ratio,BWVer, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

        xFurrow2{lauf} = x;
        yFurrow2{lauf} = y;

        Furrow2_BW_BoxStore{lauf} = BW_Box;
        Furrow2_BWHor_BoxStore{lauf} = BWHor;
        Furrow2_BWVer_BoxStore{lauf} = BWVer;
        

       % [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
       
       
        BW_Box = BW_Box_Out{3};
        BWHor = BWHor_Out{3};
        BWVer  =BWVer_Out{3};
        
        x = x_Out{3};
        y = y_Out{3};
       
      
        Furrow1Store{lauf} = Furrow1
        Furrow2Store{lauf} = Furrow2
        FurrowMeanIntensity{lauf} = (Furrow1+Furrow2)/2;

        Furrow1Raw_S1_Store(lauf) = Furrow1_S1;
        Furrow1Raw_S2_Store(lauf) = Furrow1_S2;


        Furrow2Raw_S1_Store(lauf) = Furrow2_S1;
        Furrow2Raw_S2_Store(lauf) = Furrow2_S2;

        FurrowRawMeanStore(lauf) = (Furrow2_S1 + Furrow1_S1)/2
        
         Furrow_S1_RawMeanStore(lauf) = (Furrow2_S1 + Furrow1_S1)/2
         Furrow_S2_RawMeanStore(lauf) = (Furrow2_S2 + Furrow1_S2)/2

       

        ratioFurrow = (Furrow1+Furrow2)/2;

  

        ratioFurrowStore{lauf} = ratioFurrow
       
        
        contractileRingDistance{lauf} = (pdist2([xFurrow1{lauf},yFurrow1{lauf}],[xFurrow2{lauf},yFurrow2{lauf}])) * voxelX_mum


        
end

%time = 1:timeInterval:p*timeInterval

 
anaTime = anaOnset*timeInterval;
timeMax = (p*timeInterval+p) - anaTime;
anaTime = 0;
 
 timeVec = 1:timeInterval:p*timeInterval;
            
timeVec = timeVec-timeVec(anaOnset);
    
%ingressionTime = timeVec(ingressionFrame);
             

for lauf = 1:p

    contractileRingDistanceTmp(lauf) = contractileRingDistance{lauf};
  
end

%%%%%%% Plot contractile ring diameter

h=figure
plot(timeVec,contractileRingDistanceTmp)
hold on

axis([timeVec(1) 500 0 25])  
yL = get(gca,'YLim');
%line([ingressionTime ingressionTime],yL,'Color','r');

 
xlabel ('Time [s]','FontSize', 16);
ylabel('Contractile ring dimaeter [µm]','FontSize', 16);
title(['Contractile ring diameter' tifFilename],'FontSize', 16);

print(h,'-dpdf', [tifFilename,'_FurrowIngressionDiameter.pdf']);%tifCurvetifFilename);

close all



      saveVariables = {};

            saveVariables{1} = timeVec';%time{1};
            saveVariables{2} = contractileRingDistanceTmp';
            
          
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, contractile_ring_distance'];
            
            outid = fopen([tifFilename,'ContractileRingDiameter.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'ContracitleRingDiameter.csv'],csvData,'roffset',1,'-append')

          % save([tifFilename,'.mat'])
            tifFilename
            %contractileRingDistanceTmp