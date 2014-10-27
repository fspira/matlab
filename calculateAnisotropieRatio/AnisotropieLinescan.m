%%%%% Anisotropie Linescan


clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

blue = zeros(256,3);
blue(:,3) = distvec;
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


javaaddpath '/Applications/Fiji.app/jars/ij-1.49h.jar'
addpath('/Applications/Fiji.app/scripts')
addpath('/Applications/Fiji.app/plugins/')
%addpath('/Users/spira/Desktop/programme/miji/mij.jar')
%addpath('/Users/spira/Desktop/programme/miji/ij.jar')

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/mij.jar'

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/ij.jar'


 MIJ.start('/Applications/Fiji.app')


anaOnset =1;
ingressionFrame =4;
lastFrameToConsider =11;

curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell9_100nm_vertical.lsm';
%tifFilenameMid = 'cell6_conv.tif';

saveFileName = 'cell9_100nm_vertical_FlankingCenter';

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


folderName = [tifFilename(1:truncName-1),'FlankingCenter'];
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
    
    

    imgOrigGreenTmp = img(:,:,( planeSelector));
    imgOrigGreen(:,:,lauf) = imgOrigGreenTmp{3};
    
    
    imgOrigRedTmp = img(:,:,( planeSelector));
    imgOrigRed(:,:,lauf) = imgOrigRedTmp{4};
   
     
    imgMidTmp = img(:,:,( planeSelector));
    imgMid(:,:,lauf) = imgMidTmp{1};
    
    planeSelector = planeSelector+2
    
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



%%%%% Segment the image

t3Store = doSegmentImage(S1);

[m n p] = size(t3Store)
flank1Store = {};
flank2Store = {};
flank2Mid = [];
flank1Mid = [];

%%%%%% Detect Flanks

for lauf = 1:p
    

    [flank1 flank2 flankMid1 flankMid2] = doSplitWatershedInTwo(t3Store(:,:,lauf), greenNorm(:,:,lauf))

    flank1Store{lauf} = flank1
    flank2Store{lauf} = flank2
    flank1Mid(lauf) = flankMid1
    flank2Mid(lauf) = flankMid2
    
end
    

for lauf = 1:p %length(imgGauss)
    
    flank1TmpUpdate =[];
    flank2TmpUpdate =[];
    
    
    flank1Tmp = flank1Store{lauf};
    
    flank1FirstFrame = flank1Store{1};

    flank2FirstFrame = flank2Store{1};
    
    flank1FirstCenter = flank1Mid(1);
    flank2FirstCenter = flank2Mid(1);
    
    
   
    
    flank2Tmp = flank2Store{lauf};
    
    %%%% To order the  segmented flanks 
    
    if lauf == p 
        
            flank1Tmp1 = flank1Store{lauf-1};
    
            flank2Tmp1 = flank2Store{lauf-1};
        
            flank1Dist = pdist2([ flank1FirstFrame(flank1FirstCenter(1),:)],[flank1Tmp(flank1Mid(lauf),:)])

            flank2Dist = pdist2([ flank2FirstFrame(flank2FirstCenter(1),:)],[flank2Tmp(flank2Mid(lauf),:)])

            if flank1Dist >flank2Dist

                flank1Tmp = flank2Store{lauf};
                flank2Tmp = flank1Store{lauf};

                flank1MidTmp = flank1Mid(lauf)
                flank2MidTmp = flank2Mid(lauf)

                flank1Mid(lauf) = flank2MidTmp;
                flank2Mid(lauf) = flank1MidTmp;

            else

            end  
            

    else
        
          flank1Dist = pdist2([ flank1FirstFrame(flank1FirstCenter(1),:)],[flank1Tmp(flank1Mid(lauf),:)])

            flank2Dist = pdist2([ flank1FirstFrame(flank1FirstCenter(1),:)],[flank2Tmp(flank2Mid(lauf),:)])

            if flank1Dist  > flank2Dist

                flank1Tmp = flank2Store{lauf};
                flank2Tmp = flank1Store{lauf};

                flank1MidTmp = flank1Mid(lauf)
                flank2MidTmp = flank2Mid(lauf)

                flank1Mid(lauf) = flank2MidTmp;
                flank2Mid(lauf) = flank1MidTmp;

            else

            end
            
    end
    
    MIJ.createImage(redNorm(:,:,lauf));
    MIJ.setRoi( [flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
   
    k = waitforbuttonpress 
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    flank1TmpUpdate(:,1) = tmp(:,1);
    flank1TmpUpdate(:,2) = tmp(:,2);
    
    flank1Store{lauf} = [flank1TmpUpdate(:,1) flank1TmpUpdate(:,2)]
    
    
     MIJ.setRoi( [flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
    k = waitforbuttonpress 
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    flank2TmpUpdate(:,1) = tmp(:,1);
    flank2TmpUpdate(:,2) = tmp(:,2);
    
    flank2Store{lauf} = [flank2TmpUpdate(:,1) flank2TmpUpdate(:,2)]
    
    
    MIJ.run('closeAllWindows');
    lauf

end

%%%%% Detect Poles


%for lauf = 1:p
    

 %   [pole1 pole2 poleMid1 poleMid2] = doSplitWatershedInTwoPoles(t3Store(:,:,lauf), greenNorm(:,:,lauf))

%    pole1Store{lauf} = pole1
%    pole2Store{lauf} = pole2
%    pole1Mid(lauf) = poleMid1
 %   pole2Mid(lauf) = poleMid2
    
%end

%%%%% Mark the ingressing cleavage furrow

    
           for lauf = 1:p %length(imgGauss)
     
      flank1Tmp = flank1Store{lauf};
      flank2Tmp = flank2Store{lauf};
    
    %%%%% Mark ingressing furrow
    
               

                MIJ.createImage(imgMerge(:,:,:,lauf));
                MIJ.run('8-bit');
              
               MIJ.run('Stack to RGB');
               % MIJ.selectWindow('Import from Matlab');
                MIJ.setRoi( [flank1Tmp(flank1Mid(lauf),1);flank1Tmp(flank1Mid(lauf),2)], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cNewMid1N(lauf) = coords(1);
                rNewMid1N(lauf) = coords(2);
              
             

               MIJ.setRoi( [flank2Tmp(flank2Mid(lauf),1);flank2Tmp(flank2Mid(lauf),2)], ij.gui.Roi.POINT); 
                     k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cNewMid2N(lauf) = coords(1);
                rNewMid2N(lauf) = coords(2);
           
   
    

                MIJ.run('closeAllWindows');

           end

           
           
              cNewMid1 =  cNewMid1N;
                rNewMid1 =  rNewMid1N;
                
                
                cNewMid2 =  cNewMid2N;
                rNewMid2 =  rNewMid2N;

                
                
AnalysisEnd = p;


%%%%%% Calculate cleavage furrow diameter

  for lauf = 1:  AnalysisEnd
               
            
                   dist1 = pdist2([cNewMid1(lauf),rNewMid1(lauf)],[cNewMid2(lauf),rNewMid2(lauf)]);
                   
                   distMinIngressing(lauf) = dist1*voxelX_mum;
                   
                   
                 
              
  end
            
  
    anaTime = anaOnset*timeInterval;
             timeMax = (AnalysisEnd*timeInterval+AnalysisEnd) - anaTime;
             
            anaTime = 0;
            
           
            
            
            
            timeVec = 1:timeInterval:AnalysisEnd*timeInterval;
            
            timeVec = timeVec-timeVec(anaOnset);
            
            ymax = max(distMinIngressing)+2;
            
            %%% Plot cleavage furrow diameter during ingression
            
            
           h=figure(1)
            plot(timeVec,distMinIngressing(1:AnalysisEnd)')
             axis([timeVec(1) timeMax 1 ymax])      
             set(gca,'FontSize',16,'FontName', 'Arial')

      
           
            xlabel ('Time [s]','FontSize', 16,'FontName', 'Arial');
            ylabel('Cleavage furrow diameter [µm]','FontSize', 16,'FontName', 'Arial');
            title(['Cleavage furrow ingression diameter Exp: 2433 ' tifFilename],'FontSize', 16,'FontName', 'Arial');
          %  orient landscape;
            print(h,'-dpdf', [curdir '/'  tifFilename,'_CleavageFurrowDiameter.pdf']);%tifCurvetifFilename);
            close all
            
            
            
            saveVariables = {};
 
            saveVariables{1} =timeVec';
            saveVariables{2} = distMinIngressing';
            %%%%%%%% This part save the variables
            %%%%%%%% First Variables
 
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, distance'];
            
            outid = fopen('CleavageFurrowDiameter.csv', 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ('CleavageDiameter.csv',csvData,'roffset',1,'-append')
            

    %%%%%% Get the linescan from the curves

for lauf = 1:AnalysisEnd %length(redStack)
         
         flank1Tmp = flank1Store{lauf};
         flank2Tmp = flank2Store{lauf};
         
         
         MIJ.createImage(redNorm(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRed1{lauf} = MIJ.getColumn('y');
         xRed1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         MIJ.createImage(redNorm(:,:,lauf));
         MIJ.run('setLine8');
          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRed2{lauf} = MIJ.getColumn('y');
         xRed2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
         
         
         MIJ.createImage(greenNorm(:,:,lauf));
         MIJ.run('setLine8');
          MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreen1{lauf} = MIJ.getColumn('y');
         xGreen1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         MIJ.createImage(greenNorm(:,:,lauf));
         MIJ.run('setLine8');
          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreen2{lauf} = MIJ.getColumn('y');
         xGreen2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
         
end


%test = yRedMid1

%%%%%%% Identification of midsections of the ROI's. The short script plots
%%%%%%% a high pixel value into the linescan, which can be easily detected
%%%%%%% as a peak and used to find the mid regin.

%Miji;

for lauf=1:AnalysisEnd
    %lauf =4
    redStackMid = redNorm(:,:,lauf);
    
    
         flank1Tmp = flank1Store{lauf};
         flank2Tmp = flank2Store{lauf};
         
    
    cMark1 = cNewMid1(lauf);
    rMark1 = rNewMid1(lauf);
    
    cMark2 =  cNewMid2(lauf);
    rMark2 =  rNewMid2(lauf);
     
    
         redStackMid(rMark1,cMark1) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid1{lauf} = MIJ.getColumn('y');
         xRedMid1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redMax1(lauf) = find(yRedMid1{lauf}==max(yRedMid1{lauf}))
         
         redStackMid(rMark2,cMark2) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid2{lauf} = MIJ.getColumn('y');
         xRedMid2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redMax2(lauf) = find(yRedMid2{lauf}==max(yRedMid2{lauf}))
         
         
       

end


%%%%% Calculate ratio between green and red

for lauf = 1:AnalysisEnd
    
    ratioFurrow1{lauf} = yGreen1{lauf}./yRed1{lauf};
    ratioFurrow2{lauf} = yGreen2{lauf}./yRed2{lauf};


end


            

[greenOut redOut midOut] = doAverageFlanks(redMax1,redMax2,yGreen1,yGreen2,yRed1,yRed2)


for lauf = 1:AnalysisEnd
    
    ratioFurrow1Merged{lauf} = greenOut{lauf}./redOut{lauf};
    


end

%%%%%% The following section shifts the arrays of the linescan in a way
%%%%%% that every linescan has 300px in lenght and is centered.
for lauf = 1: length(midOut)
    if midOut(lauf) < 150


        linescanSelectorTmp_I =  ratioFurrow1Merged{lauf}

       % length(linescanSelectorTmp_I);
       % midOut(lauf);

        addZero = zeros(1,150 - midOut(lauf))';

        linescanSelectorCrop = cat(1,addZero,linescanSelectorTmp_I);
        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorCrop,addZero);

        linescanSelectorCropStoreGreen(:,lauf) = linescanSelectorCrop(1:300);
    elseif midOut(lauf) > 150

         linescanSelectorTmp_I = ratioFurrow1Merged{lauf}


        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreen(:,lauf) =  linescanSelectorCrop(midOut(lauf)-149:300+(midOut(lauf)-150));


    elseif midOut(lauf) == 150

        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1, linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreen(:,lauf) = linescanSelectorCrop(1:300);

    end
end

%%%%%% The following section shifts the arrays of the linescan in a way
%%%%%% that every linescan has 300px in lenght and is centered.

longDistanceVec = ((-149:150)*voxelX_mum)'

imshow(linescanSelectorCropStoreGreen,[])


colormap(jet)

normalizedKymo = normalizedImage(linescanSelectorCropStoreGreen);

imwrite(normalizedKymo,jet,[tifFilename,'Aniso_Kymo_colormapJet.tif']);
imwrite(normalizedKymo,[tifFilename,'Aniso_Kymo.tif']);

imwrite(A,map,filename)

 clear linescanSelectorTmp
    linescanSelectorTmp = zeros(300,43);
     linescanSelectorTmp  = double(linescanSelectorTmp)
midOut


linescanSelectorIn = redOut;

[greenOverlayStore]= doGenerateAnisotropieOverlay(greenNorm, ratioFurrow1, ratioFurrow2,flank1Store,flank2Store,AnalysisEnd);
  
         tiffwrite_mat(greenOverlayStore,[tifFilename,'Aniso_Overlay.tif']);
         
       greenOverlayStoreNorm =   normalizedImage3D(greenOverlayStore);
                 
       tiffwrite_mat(greenOverlayStoreNorm,[tifFilename,'Aniso_OverlayNorm.tif']);
       

       for lauf = 1:AnalysisEnd
        
       plot(1:length(linescanSelectorCropStoreGreen(:,lauf)),linescanSelectorCropStoreGreen(:,lauf))
       hold on
       
       end

    %%%%%%%%% Normalize data by subtracting 5% of flanking pixels first
    %%%%%%%%% selecting a 15µm window around the central pixel
    
   
    
    windowBoundaries = round(10/voxelX_mum);
    
for lauf = 1:AnalysisEnd

   linescanSelectorToCrop = linescanSelectorCropStoreGreen(:,lauf);
    linescanSelectorTmp = linescanSelectorToCrop(150-(windowBoundaries-1):(150+windowBoundaries))
   
    
    averagePixel = round((length(linescanSelectorTmp)/100)*5);
    linescanAvgLeft = linescanSelectorTmp(1:averagePixel);
    linescanAvgRight =linescanSelectorTmp(end-averagePixel:end);
    
    linescanSubMean = mean((linescanAvgLeft)+mean(linescanAvgRight))/2;
    
 linescanSelector{lauf} = linescanSelectorTmp - linescanSubMean; 
  
 %linescanSelector{lauf} =  linescanSelectorTmp
   
end

    newDistanceVec = 150-(windowBoundaries-1):(150+windowBoundaries)
    newDistance = length(newDistanceVec)/2;
    
    timeVecTmp = timeInterval:timeInterval:timeInterval*AnalysisEnd;
    
   % timeVec = timeVecTmp - (timeInterval * anaOnset)
   timeVec  =timeVecTmp;

%%%% Use the average of the total intensity for normalization
%%%% Has to be implemented
for lauf = 1:length(linescanSelectorIn)
 

    
    %%%% find frame relative to anaphase onset
   % time1 = timeVec(lauf);
    
    
    
    %%%% load the data from the linescan and find the center position
    redMaxTmp = (0-(newDistance-1):newDistance)*voxelX_mum;
    yGreenTmp = linescanSelector{lauf}%- min(linescanSelector{lauf})
    %%%% Calculated the full width half maximum
    
    %xValue = inputLinescan(zeroPos-10:zeroPos+10);
   
    
    maxX =yGreenTmp(newDistance)
  
    newY =newDistance;
    
    ft = fittype( 'gauss1' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.Robust = 'Bisquare';
   % opts.StartPoint = [maxX newY 0.5];
    opts.Upper = [Inf Inf Inf];
    
     [fittedmodel, goodness, output] = fit( redMaxTmp',yGreenTmp,ft,opts);
     
      sigma1Orig = fittedmodel.c1/sqrt(2); %%%% Gives sigma
          FWHMOrig = 2*sqrt(2*log(2))*sigma1Orig; %%% Calculates the FWHM

        x_axOrig = -400:0.01:400; % axis
        y_fittedOrig = fittedmodel(x_axOrig); 
       

    
   
    FWHMOrigStore(lauf) = FWHMOrig;

 %   h=figure(1)
    %%%% Plot the data
  %  plot(((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum),yGreenTmp,'c')
    %X = [((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum)',yGreenTmp]
    %legend([num2str(time1) ' s'])
   % hold on
       
    %   close all
    x_axOrig_store{lauf} = x_axOrig;
    y_fittedOrig_store{lauf} = y_fittedOrig;
  %  plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
end
     FWHMOrigStore =  FWHMOrigStore'

curdir = pwd;
    mkdir('redfit')
    
    for lauf = 1:AnalysisEnd
              h = figure(1) 
        plot((0-(newDistance-1):newDistance)*voxelX_mum,linescanSelector{lauf}')
       
     
      
        hold on
          set(gca,'FontSize',16,'FontName', 'Arial')
      plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
       axis([-15 +15 -1 1]) 
  title('Anisotropie SiR-actin green/red Exp. 2433','FontSize', 16,'FontName', 'Arial')
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial');
    ylabel('Ratio green/red  [A.U.]' ,'FontSize', 16,'FontName', 'Arial');
    legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(timeVec(lauf)) 's' ])
        pause(0.2)
       
       
      % print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5percAVG.pdf']);%tifCurvetifFilename);
   
      print(h,'-dpdf', [curdir '\' ,'redFit' '\' , 'frame_', num2str(lauf) ,'_GaussFit_FWHM_Anisotropie.pdf']);%tifCurvetifFilename);
   close all
    end
    
    
    
         
    for lauf = 1:AnalysisEnd
              h = figure(1) 
        plot(((0-149):150)*voxelX_mum, linescanSelectorCropStoreGreen(:,lauf)')
       
     
      
        hold on
          set(gca,'FontSize',16,'FontName', 'Arial')
     % plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
       axis([-25 +25 0 1.5]) 
  title('Anisotropie SiR-actin green/red Exp. 2433','FontSize', 16,'FontName', 'Arial')
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial');
    ylabel('Ratio green/red  [A.U.]' ,'FontSize', 16,'FontName', 'Arial');
    legend( [num2str(timeVec(lauf)), 's'] )
        pause(0.2)
       
       
      % print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5percAVG.pdf']);%tifCurvetifFilename);
   
      print(h,'-dpdf', [curdir '\' ,'redFit' '\' , 'frame_', num2str(lauf) ,'_FullLinescanLength_Anisotropie.pdf']);%tifCurvetifFilename);
   close all
    end
    