
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


anaOnset = 1;
ingressionFrame = 4;
lastFrameToConsider = 8;

curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell3_100nm_horizontal.lsm';
%tifFilenameMid = 'cell6_conv.tif';

saveFileName = 'cell3_100nm_horizontal_center';

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


for lauf = 1:p
        
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
       
       xPole1{lauf} = x;
       yPole{lauf} = y;
       
       [Pole1, Pole1_S1, Pole1_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box, S1(:,:,lauf), S2(:,:,lauf))

        

        Pole1_BW_BoxStore{lauf} = BW_Box;
        Pole1_BWHor_BoxStore{lauf} = BWHor;
        Pole1_BWVer_BoxStore{lauf} = BWVer;
        
        
         BW_Box = BW_Box_Out{4};
        BWHor = BWHor_Out{4};
        BWVer  =BWVer_Out{4};
        
        x = x_Out{4};
        y = y_Out{4};

    %    [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge(:,:,:,lauf),analysisFrame,ratio,voxelX_mum);
        [Pole2, Pole2_S1, Pole2_S2] = doCalculateRatioBoundingBox(ratio,BWHor, BW_Box,S1(:,:,lauf), S2(:,:,lauf))

        xPole2{lauf} = x;
        yPole2{lauf} = y;

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
        
         Furrow_S1_RawMeanStore(lauf) = (Furrow2_S1 + Furrow1_S1)/2
         Furrow_S2_RawMeanStore(lauf) = (Furrow2_S2 + Furrow1_S2)/2

        Pole1Raw_S1_Store(lauf) = Pole1_S1
        Pole1Raw_S2_Store(lauf) = Pole1_S2

        Pole2Raw_S1_Store(lauf) = Pole2_S1
        Pole2Raw_S2_Store(lauf) = Pole2_S2

        PoleRawMeanStore(lauf) = (Pole1_S1 + Pole1_S2)/2
         
        Pole_S1_RawMeanStore(lauf) = (Pole1_S1 + Pole1_S2)/2
        Pole_S2_RawMeanStore(lauf) = (Pole1_S2 + Pole1_22)/2


        
        Pole1Store{lauf} = Pole1
        Pole2Store{lauf} = Pole2
        PoleMeanIntensity{lauf}= (Pole1+Pole2)/2;

        ratioFurrow = (Furrow1+Furrow2)/2;

        ratioPole = (Pole1+Pole2)/2;

        ratioFurrowStore{lauf} = ratioFurrow
        ratioPoleStore{lauf} = ratioPole
        
        
        contractileRingDistance{lauf} = (pdist2([xFurrow1{lauf},yFurrow1{lauf}],[xFurrow2{lauf},yFurrow2{lauf}])) * voxelX_mum


        
end

%time = 1:timeInterval:p*timeInterval

 
anaTime = anaOnset*timeInterval;
timeMax = (p*timeInterval+p) - anaTime;
anaTime = 0;
 
 timeVec = 1:timeInterval:p*timeInterval;
            
timeVec = timeVec-timeVec(anaOnset);
    
ingressionTime = timeVec(ingressionFrame);
             

for lauf = 1:p

    contractileRingDistanceTmp(lauf) = contractileRingDistance{lauf};
    ratioFurrowStoreTmp(lauf) =  ratioFurrowStore{lauf}
    ratioPoleStoreTmp(lauf) = ratioPoleStore{lauf}
    FurrowMeanIntensityTmp(lauf) = FurrowMeanIntensity{lauf}
    PoleMeanIntensityTmp(lauf) = PoleMeanIntensity{lauf}

end

%%%%%%% Plot contractile ring diameter

h=figure
plot(timeVec,contractileRingDistanceTmp)
hold on

axis([timeVec(1) 500 0 25])  
yL = get(gca,'YLim');
line([ingressionTime ingressionTime],yL,'Color','r');

 
xlabel ('Time [s]','FontSize', 16);
ylabel('Contractile ring dimaeter [µm]','FontSize', 16);
title(['Contractile ring diameter' tifFilename],'FontSize', 16);

print(h,'-dpdf', [tifFilename,'_FurrowIngressionDiameter.pdf']);%tifCurvetifFilename);

close all

%%%%%% plot ratio intensities

h=figure

 plot(timeVec, FurrowMeanIntensityTmp,'-r')
 hold on
  plot(timeVec, PoleMeanIntensityTmp,'-b')

axis([timeVec(1) 500 0.7 1.6])   

yL = get(gca,'YLim');
line([ingressionTime ingressionTime],yL,'Color','r');


legend('ratio furrow', ...
  'ratio pole')


xlabel ('Time [s]','FontSize', 16);
ylabel('Ratio [A.U.]','FontSize', 16);
title(['Ratio green/red furrow and pole' tifFilename],'FontSize', 16);

print(h,'-dpdf', [tifFilename,'_ratioFurrowPole.pdf']);%tifCurvetifFilename);


close all

%%%%%% Plot raw intensities


h = figure


 plot(timeVec, FurrowRawMeanStore,'-r')
 hold on
  plot(timeVec,PoleRawMeanStore,'-b')

axis([timeVec(1) 500 0 1400])   

  
yL = get(gca,'YLim');
line([ingressionTime ingressionTime],yL,'Color','r');


legend('raw furrow', ...
  'raw pole')


xlabel ('Time [s]','FontSize', 16);
ylabel('Raw inensities [A.U.]','FontSize', 16);
title(['Intensities furrow and pole' tifFilename],'FontSize', 16);

print(h,'-dpdf', [tifFilename,'_rawFurrowPole.pdf']);%tifCurvetifFilename);


close all



% 
%     
% x1 = timeVec;
% y1 = contractileRingDistance;
% y2 = FurrowMeanIntensityTmp;
% y3 = PoleMeanIntensityTmp;
% y4 = FurrowRawMeanStore;
% y5 = PoleRawMeanStore;
% 
% h= figure
% [AX,H1,H2] = plotyy(x1,y2,x1,y3,x1,y4, 'plot');
%   
%   set(get(AX(1),'Ylabel'),'String','Distance [µm]','FontSize', 16) 
%  %  set(get(AX(1),'Ylabel'),'String','Distance [µm]','FontSize', 16) 
%   set(get(AX(2),'Ylabel'),'String','Normalized mean intensity [A.U.]','FontSize', 16) 
%        
% hold(AX(1), 'on')
% hold(AX(2), 'on')
% %timeActinin, ActininAna
% % Plot the third curve
% h3 = plot(timeMRLC,LifeactAna,'r', 'Parent', AX(2));
% h4 = plot(timeMRLC,MRLCAna,'c', 'Parent', AX(2));   
% h5 = plot(timeRhoA,RhoAAna,'b', 'Parent', AX(2));
% h5 = plot(timeActinin(1:36), ActininAna(1:36),'g', 'Parent', AX(2)); 
% 
% 
%             yL = get(gca,'YLim');
%             line([120 120],yL,'Color','r');
%             line([250 250],yL,'Color','r');

  
  
  

      saveVariables = {};

            saveVariables{1} = timeVec';%time{1};
            saveVariables{2} = contractileRingDistanceTmp';
            
            
            saveVariables{3} = ratioFurrowStoreTmp';
            saveVariables{4} = ratioPoleStoreTmp';
            saveVariables{5} = FurrowMeanIntensityTmp'
            saveVariables{6} = PoleMeanIntensityTmp';
            saveVariables{7} = FurrowRawMeanStore';
            saveVariables{8} = PoleRawMeanStore';
            
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, contractile_ring_distance, ratio_furrow, ratio_pole',...
                'furrow_mean_intensity, pole_mean_intensity, furrow raw, pole raw'];
            
            outid = fopen([tifFilename,'RatioAnalysis_mid.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([tifFilename,'RatioAnalysis_mid.csv'],csvData,'roffset',1,'-append')

           save([tifFilename,'.mat'])
            tifFilename