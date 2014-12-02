%%%%% anisotropy Linescan


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
addpath('/Users/spira/Desktop/programme/calculateanisotropyRatio')
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


anaOnset =14;
ingressionFrame =19;
lastFrameToConsider =31;
AnalysisEnd = lastFrameToConsider;

curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell2_100nm_vertical.lsm';
%tifFilenameMid = 'cell6_conv.tif';

saveFileName = 'cell2_100nm_vertical_entireFlank';

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


greenImg = uint16(imgOrigGreen);
redImg = uint16(imgOrigRed);


[greenImg redImg] = bkgCorrectionRedGreen(greenImg, redImg);

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

%%%%% Calculate anisotropy

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
    
%%%%%% Burn ROI's into the image
flank1Tmp = flank1Store{5};
flank2Tmp = flank2Store{5};
flank1Tmp = pole1Store{5};
flank1Tmp = pole2Store{5};

templateImage = S1(:,:,5);
templateImage(:,:) = 0;
    
for lauf = 1:length(flank1Tmp)
  
    templateImage(flank1Tmp(lauf,2),flank1Tmp(lauf,1)) = 255;
    
end



imwrite(templateImage,'splitSegment_cell9_Frame5_Pole2.tif');
counter = 26
for lauf = counter:AnalysisEnd
                
     flank1Store{lauf} = [50 50]; 
    flank2Store{lauf} = [50 50];
    flank1Mid(lauf) = 1
    flank2Mid(lauf) = 1
    
end

for lauf =1:AnalysisEnd
    MIJ.closeAllWindows
    flank1TmpUpdate =[];
    flank2TmpUpdate =[];
    
    
    flank1Tmp = flank1Store{lauf};
    
    flank1FirstFrame = flank1Store{1};

    flank2FirstFrame = flank2Store{1};
    
    flank1FirstCenter = flank1Mid(1);
    flank2FirstCenter = flank2Mid(1);
    
    
   
    
    flank2Tmp = flank2Store{lauf};
    
    %%%% To order the  segmented flanks 
    
    if lauf == AnalysisEnd
        
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
    clear logData
     MIJ.run('getSplineCoordinates')
     logData = MIJ.getLog;
     MIJ.run('closeLogWindow')
    logDataTmp = char(logData);
    logDataTmp1 = strread(logDataTmp, '%s');
        
         xLogTmp = logDataTmp1(2:4:end);
         yLogTmp = logDataTmp1(3:4:end);
         
             
        xLog = zeros(length(xLogTmp),1,'double');
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
          flank1StoreNew{lauf} = [xLog yLog]
    
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    flank1TmpUpdate(:,1) = xLog;
    flank1TmpUpdate(:,2) = yLog;
    
    flank1Store{lauf} = [flank1TmpUpdate(:,1) flank1TmpUpdate(:,2)]
    
        
     MIJ.createImage(redNorm(:,:,lauf));
     MIJ.setRoi( [flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
    k = waitforbuttonpress 
     clear logData
    
     MIJ.run('getSplineCoordinates')
     logData = MIJ.getLog;
        MIJ.run('closeLogWindow')
    logDataTmp = char(logData);
    logDataTmp1 = strread(logDataTmp, '%s');
        
         xLogTmp = logDataTmp1(2:4:end);
         yLogTmp = logDataTmp1(3:4:end);
         
             
        xLog = zeros(length(xLogTmp),1,'double');
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
          flank2StoreNew{lauf} = [xLog yLog]
        
    
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    flank2TmpUpdate(:,1) = xLog;
    flank2TmpUpdate(:,2) = yLog;
    
    
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
  %  flank2TmpUpdate(:,1) = tmp(:,1);
  %  flank2TmpUpdate(:,2) = tmp(:,2);
      MIJ.closeAllWindows
    flank2Store{lauf} = [flank2TmpUpdate(:,1) flank2TmpUpdate(:,2)]
    
    
    MIJ.run('closeAllWindows');
    lauf

end






%%%%% Detect Poles

counter = 1;
for lauf = 1:AnalysisEnd
    

    [pole1 pole2 poleMid1 poleMid2] = doSplitWatershedInTwoPoles(t3Store(:,:,lauf), greenNorm(:,:,lauf))

    pole1Store{lauf} = pole1
    pole2Store{lauf} = pole2
    pole1Mid(lauf) = poleMid1
    pole2Mid(lauf) = poleMid2
    counter = counter+1;
end

for lauf = counter:AnalysisEnd
    
    
    pole1Store{lauf} = [50 50]; 
    pole2Store{lauf} = [50 50]; 
     pole1Mid(lauf) = 1
     pole2Mid(lauf) =  1

end

%%%%% This part orderes the segmented poles
for lauf = 1:AnalysisEnd %length(imgGauss)
    
    pole1TmpUpdate =[];
    pole2TmpUpdate =[];
    
    
    pole1Tmp = pole1Store{lauf};
    
    pole1FirstFrame = pole1Store{1};

    pole2FirstFrame = pole2Store{1};
    
    pole1FirstCenter = pole1Mid(1);
    pole2FirstCenter = pole2Mid(1);
    
    
   
    
    pole2Tmp = pole2Store{lauf};
    
    %%%% To order the  segmented poles 
    
    if lauf == AnalysisEnd
        
            pole1Tmp1 = pole1Store{lauf-1};
    
            pole2Tmp1 = pole2Store{lauf-1};
        
            pole1Dist = pdist2([ pole1FirstFrame(pole1FirstCenter(1),:)],[pole1Tmp(pole1Mid(lauf),:)])

            pole2Dist = pdist2([ pole2FirstFrame(pole2FirstCenter(1),:)],[pole2Tmp(pole2Mid(lauf),:)])

            if pole1Dist >pole2Dist

                pole1Tmp = pole2Store{lauf};
                pole2Tmp = pole1Store{lauf};

                pole1MidTmp = pole1Mid(lauf)
                pole2MidTmp = pole2Mid(lauf)

                pole1Mid(lauf) = pole2MidTmp;
                pole2Mid(lauf) = pole1MidTmp;

            else

            end  
            

    else
        
          pole1Dist = pdist2([ pole1FirstFrame(pole1FirstCenter(1),:)],[pole1Tmp(pole1Mid(lauf),:)])

            pole2Dist = pdist2([ pole1FirstFrame(pole1FirstCenter(1),:)],[pole2Tmp(pole2Mid(lauf),:)])

            if pole1Dist  > pole2Dist

                pole1Tmp = pole2Store{lauf};
                pole2Tmp = pole1Store{lauf};

                pole1MidTmp = pole1Mid(lauf)
                pole2MidTmp = pole2Mid(lauf)

                pole1Mid(lauf) = pole2MidTmp;
                pole2Mid(lauf) = pole1MidTmp;

            else

            end
            
    end
    
    MIJ.createImage(redNorm(:,:,lauf));
    MIJ.setRoi( [pole1Tmp(:,1)'; pole1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
   
    k = waitforbuttonpress 
    
     clear logData
     MIJ.run('getSplineCoordinates')
     logData = MIJ.getLog;
     MIJ.run('closeLogWindow')
    logDataTmp = char(logData);
    logDataTmp1 = strread(logDataTmp, '%s');
        
         xLogTmp = logDataTmp1(2:4:end);
         yLogTmp = logDataTmp1(3:4:end);
         
             
        xLog = zeros(length(xLogTmp),1,'double');
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
         Pole1StoreNew{lauf} = [xLog yLog]
    
    
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    pole1TmpUpdate(:,1) = xLog;
    pole1TmpUpdate(:,2) = yLog;
    
    pole1Store{lauf} = [pole1TmpUpdate(:,1) pole1TmpUpdate(:,2)]
    
     MIJ.createImage(redNorm(:,:,lauf));
     MIJ.setRoi( [pole2Tmp(:,1)'; pole2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
    k = waitforbuttonpress 
    
    clear logData
     MIJ.run('getSplineCoordinates')
     logData = MIJ.getLog;
     MIJ.run('closeLogWindow')
    logDataTmp = char(logData);
    logDataTmp1 = strread(logDataTmp, '%s');
        
         xLogTmp = logDataTmp1(2:4:end);
         yLogTmp = logDataTmp1(3:4:end);
         
             
        xLog = zeros(length(xLogTmp),1,'double');
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
         Pole2StoreNew{lauf} = [xLog yLog]
    
    
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    pole2TmpUpdate(:,1) = xLog;
    pole2TmpUpdate(:,2) = yLog;
    
    
    pole2Store{lauf} = [pole2TmpUpdate(:,1) pole2TmpUpdate(:,2)]
    
    
    MIJ.run('closeAllWindows');
    lauf

end



%

%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%% Mark the ingressing cleavage furrow

    
           for lauf = 1:AnalysisEnd %length(imgGauss)
     
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

                
                for lauf = 1:AnalysisEnd
                
                    imshow(imgMerge(:,:,:,lauf))
                    hold on
                    plot(cNewMid1(lauf),rNewMid1(lauf),'xb')
                    plot(cNewMid2(lauf),rNewMid2(lauf),'xr')
                    pause(0.5)
                    
                end
                
         



%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%% Mark the poles

    
           for lauf = 1:AnalysisEnd %length(imgGauss)
     
      flank1Tmp = pole1Store{lauf};
      flank2Tmp = pole2Store{lauf};
    
    %%%%% Mark ingressing furrow
    
               

                MIJ.createImage(imgMerge(:,:,:,lauf));
                MIJ.run('8-bit');
              
               MIJ.run('Stack to RGB');
               % MIJ.selectWindow('Import from Matlab');
                MIJ.setRoi( [flank1Tmp(pole1Mid(lauf),1);flank1Tmp(pole1Mid(lauf),2)], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cPoleMid1N(lauf) = coords(1);
                rPoleMid1N(lauf) = coords(2);
              
             

               MIJ.setRoi( [flank2Tmp(pole2Mid(lauf),1);flank2Tmp(pole2Mid(lauf),2)], ij.gui.Roi.POINT); 
                     k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cPoleMid2N(lauf) = coords(1);
                rPoleMid2N(lauf) = coords(2);
           
   
    

                MIJ.run('closeAllWindows');

           end

lauf = 5
          templateImage = normalizedImage(S1(:,:,5));
          
          templateImage(rPoleMid1N(lauf),cPoleMid1N(lauf)) = 255;
          templateImage(rPoleMid2N(lauf),cPoleMid2N(lauf)) = 255;
          templateImage( rNewMid1(lauf), cNewMid1(lauf)) = 255;
          templateImage( rNewMid2(lauf), cNewMid2(lauf)) = 255;
          
          imshow(templateImage,[])
           


imwrite(templateImage,'midPoints_cell9_Frame5_Pole2.tif');



imwrite(BW2,'segmented_Chromatin_cell2_Frame5.tif');

[cChromatin1 rChromatin1 cChromatin2 rChromatin2] = doChromatinChromatinDistance(imgMid)

 chromatin1Mid =  zeros(AnalysisEnd,2);
 
 chromatin2Mid = zeros(AnalysisEnd,2);

for lauf = 1:AnalysisEnd
    
     chromatin1Mid(lauf,1:2) = [cChromatin1(lauf) rChromatin1(lauf)]
     chromatin2Mid(lauf,1:2) = [cChromatin2(lauf) rChromatin2(lauf)]
    

end

[chromatin1Mid chromatin2Mid] = doSortChromatinPositions(chromatin1Mid, chromatin2Mid, AnalysisEnd)


cChromatin1 = chromatin1Mid(:,1);
rChromatin1 = chromatin1Mid(:,2);

cChromatin2 = chromatin2Mid(:,1);
rChromatin2 = chromatin2Mid(:,2);


            %%%%
            
            for lauf = 1:AnalysisEnd
                
                MIJ.createImage(imgMid(:,:,lauf));
                MIJ.run('8-bit');
              
               % MIJ.selectWindow('Import from Matlab');
                MIJ.setRoi( [cChromatin1(lauf);rChromatin1(lauf)], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cChromatin1(lauf) = coords(1);
                rChromatin1(lauf) = coords(2);
            
                MIJ.setRoi( [cChromatin2(lauf);rChromatin2(lauf)], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
               cChromatin2(lauf) = coords(1);
               rChromatin2(lauf) = coords(2);
               
               


                MIJ.run('closeAllWindows');


               
            end

%%%%%% Calculate cleavage furrow diameter

  for lauf = 1:  AnalysisEnd
               
            
                   dist1 = pdist2([cNewMid1(lauf),rNewMid1(lauf)],[cNewMid2(lauf),rNewMid2(lauf)]);
                   
                   distMinIngressing(lauf) = dist1*voxelX_mum;
                   
                   
                 
              
  end
  
  
  
%%%%%% Calculate chromatin chromatin diameter

  for lauf = 1:  AnalysisEnd
               
            
                   dist1 = pdist2([cChromatin1(lauf),cChromatin2(lauf)],[cChromatin2(lauf),rChromatin2(lauf)]);
                   
                   distChromatinChromatin(lauf) = dist1*voxelX_mum;
                   
                   
                 
              
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
            title(['Cleavage furrow diameter Exp: 2433 ' tifFilename],'FontSize', 16,'FontName', 'Arial');
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
            

    %%%%%% Generate linescan using the positions identified for flank1/2
    %%%%%% and pole1/2

    clear yRed1 xRed1 yRed2 xRed2 yGreen1 xGreen1 yGreen2 xGreen2
    
for lauf = 1:AnalysisEnd %length(redStack)
    
    %%%%%% measure linescans for flanks
    
         flank1Tmp = flank1Store{lauf};
         flank2Tmp = flank2Store{lauf};
         
         
         MIJ.createImage(redNorm(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         
         
         %%%%% Update InterpolatedLine
          clear logData
          MIJ.run('getSplineCoordinates')
          logData = MIJ.getLog;
          
          MIJ.run('closeLogWindow')
         
          MIJ.run('getLinescanRed');
         yRed1{lauf} = MIJ.getColumn('y');
         xRed1{lauf} = MIJ.getColumn('x');
        logDataTmp = char(logData);
        logDataTmp1 = strread(logDataTmp, '%s');
        
         xLogTmp = logDataTmp1(2:4:end);
         yLogTmp = logDataTmp1(3:4:end);
         
             
        xLog = zeros(length(xLogTmp),1,'double');
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
         flank1StoreNew{lauf} = [xLog yLog]
       
         
         MIJ.createImage(redNorm(:,:,lauf));
         MIJ.run('setLine8');
          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         
        
         
         %%%% Update interpolatedLine
              clear logData
              MIJ.run('getSplineCoordinates')
              logData = MIJ.getLog;

             MIJ.run('closeLogWindow')
          
            MIJ.run('getLinescanRed');
            yRed2{lauf} = MIJ.getColumn('y');
            xRed2{lauf} = MIJ.getColumn('x');
            
            
            logDataTmp = char(logData);
            logDataTmp1 = strread(logDataTmp, '%s');

             xLogTmp = logDataTmp1(2:4:end);
             yLogTmp = logDataTmp1(3:4:end);


            xLog = zeros(length(xLogTmp),1,'double');
            yLog = zeros(length(yLogTmp),1,'double');

             for subrun = 1 : length(xLogTmp)

                 xLog(subrun) = round(str2num(xLogTmp{subrun}));
                 yLog(subrun) = round(str2num(yLogTmp{subrun}));

             end

             flank2StoreNew{lauf} = [xLog yLog]
            
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
         
         
         
         
         
         
            %%%%%% Linescan at the poles
%           
%         pole1Tmp = pole1Store{lauf};
%         pole2Tmp = pole2Store{lauf};
%          
%          
%         MIJ.createImage(redNorm(:,:,lauf));
%         MIJ.run('setLine8');
%         MIJ.setRoi( [ pole1Tmp(:,1)'; pole1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%         MIJ.run('getLinescanRed');
%         yRedPole1{lauf} = MIJ.getColumn('y');
%         xRedPole1{lauf} = MIJ.getColumn('x');
%         MIJ.run('closeAllWindows');
%          
%         MIJ.createImage(redNorm(:,:,lauf));
%         MIJ.run('setLine8');
%          MIJ.setRoi( [ pole2Tmp(:,1)'; pole2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%         MIJ.run('getLinescanRed');
%         yRedPole2{lauf} = MIJ.getColumn('y');
%         xRedPole2{lauf} = MIJ.getColumn('x');
%         MIJ.run('closeAllWindows');
%          
%          
%          
%          
%         MIJ.createImage(greenNorm(:,:,lauf));
%         MIJ.run('setLine8');
%          MIJ.setRoi( [ pole1Tmp(:,1)'; pole1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%         MIJ.run('getLinescanRed');
%         yGreenPole1{lauf} = MIJ.getColumn('y');
%         xGreenPole1{lauf} = MIJ.getColumn('x');
%         MIJ.run('closeAllWindows');
%          
%         MIJ.createImage(greenNorm(:,:,lauf));
%         MIJ.run('setLine8');
%          MIJ.setRoi( [ pole2Tmp(:,1)'; pole2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%         MIJ.run('getLinescanRed');
%         yGreenPole2{lauf} = MIJ.getColumn('y');
%         xGreenPole2{lauf} = MIJ.getColumn('x');
%         MIJ.run('closeAllWindows');
         
        
         
end

flank1Old = flank1Store;
flank2Old = flank2Store;

flank1Store = flank1StoreNew;
flank2Store = flank2StoreNew;
%test = yRedMid1

%%%%%%% Identification of midsections of the ROI's. The following piece of
%%%%%%% code plots high intensity pixels at the center position into the
%%%%%%% linescan, which then can easily be detected by matlab.

%Miji;

for lauf=1:AnalysisEnd
    %lauf =4
    redStackMid = redNorm(:,:,lauf);
    redStackMid(:,:) = 0;
    
    
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
         
         redStackMid(rMark2,cMark2) = 650000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid2{lauf} = MIJ.getColumn('y');
         xRedMid2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redMax2(lauf) = find(yRedMid2{lauf}==max(yRedMid2{lauf}))
         
         
         %%%%%% Detection of the polar region
         
%          
%         pole1Tmp = pole1Store{lauf};
%         pole2Tmp = pole2Store{lauf};
%          
%     
%     cMark1 = cPoleMid1N(lauf);
%     rMark1 = rPoleMid1N(lauf);
%     
%     cMark2 =  cPoleMid2N(lauf);
%     rMark2 =  rPoleMid2N(lauf);
%        redStackMid(:,:) = 0;
%     
%          redStackMid(rMark1,cMark1) = 65000;
%         MIJ.createImage(redStackMid(:,:));
%         MIJ.run('setLine12');
%         MIJ.setRoi( [ pole1Tmp(:,1)'; pole1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%         MIJ.run('getLinescanRed');
%         yRedMidPole1{lauf} = MIJ.getColumn('y');
%         xRedMidPole1{lauf} = MIJ.getColumn('x');
%         MIJ.run('closeAllWindows');
%          
%         redMaxPole1(lauf) = find(yRedMidPole1{lauf}==max(yRedMidPole1{lauf}))
%          
%         redStackMid(rMark2,cMark2) = 65000;
%         MIJ.createImage(redStackMid(:,:));
%         MIJ.run('setLine12');
%         MIJ.setRoi( [ pole2Tmp(:,1)'; pole2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%         MIJ.run('getLinescanRed');
%         yRedMidPole2{lauf} = MIJ.getColumn('y');
%         xRedMidPole2{lauf} = MIJ.getColumn('x');
%         MIJ.run('closeAllWindows');
%          
%        redMaxPole2(lauf) = find(yRedMidPole2{lauf}==max(yRedMidPole2{lauf}))
%          
%          
         
       

end


%%%%% Calculate ratio between green and red

for lauf = 1:AnalysisEnd
    
    ratioFurrow1{lauf} = yGreen1{lauf}./yRed1{lauf};
    ratioFurrow2{lauf} = yGreen2{lauf}./yRed2{lauf};
    %ratioPole1{lauf} = yGreenPole1{lauf}./ yRedPole1{lauf};
    %ratioPole2{lauf} = yGreenPole2{lauf}./ yRedPole2{lauf};



end


            

[greenOut redOut midOut] = doAverageFlanks(redMax1,redMax2,yGreen1,yGreen2,yRed1,yRed2)


[greenOutPole redOutPole midOutPole] = doAverageFlanks(redMaxPole1,redMaxPole2,yGreenPole1,yGreenPole2,yRedPole1,yRedPole2)


%%%% Normalize between 0 and 1

for lauf = 1:AnalysisEnd
    
    ratioFurrow1Merged{lauf} = greenOut{lauf}./redOut{lauf};
   %ratioPoleMerged{lauf} = greenOutPole{lauf}./redOutPole{lauf};
   


end


    ratioTmp = ratioFurrow1Merged{16}
    ratioFurrowMean = mean(ratioTmp(redMax1(16)-10:redMax1(16)+10))
  
    
   
 %   ratioPoleMean = mean(ratioTmp(220-10:220+10))
    
    redMaxPole2 = 125;
    
  ratioTmpPole = ratioPoleMerged{1}
  ratioPoleMean = mean(ratioTmp(260-10:260+10))

    
    %%%%% Normalize images between min and max - relative to the first
    %%%%% frame. Min intensities are seen at the ingressing cleavage furrow
    %%%%% and max intensities at the poles

for lauf =1 :AnalysisEnd
    
    ratioTmp = ratioFurrow1Merged{lauf}
    ratioNorm = ((ratioTmp- ratioFurrowMean)/(ratioPoleMean- ratioFurrowMean)) 
 %    ratioTmpPole = ratioPoleMerged{lauf}
  %  ratioNormPole = ((ratioTmpPole- ratioFurrowMean)/(ratioPoleMean- ratioFurrowMean)) 
    
    
end


%%%%%% Caculations for the furrow
    ratioTmp1 = ratioFurrow1{16}
    ratioTmp2 = ratioFurrow1{16}
    
  %  ratioTmpPole1 = ratioPole1{1}
  %  ratioTmpPole2 = ratioPole2{1}
   
    ratioFurrowMean1 = min(ratioTmp1)
    ratioFurrowMean2= min(ratioTmp2)
    ratioPoleMean1 = max(ratioTmp1)
    ratioPoleMean2 = max(ratioTmp1)%%%%% has to be changed to pole for the final version
   
  
    
    
ratioNorm = {}
ratioNormFlank1 = {}
ratioNormFlank2 = {}
for lauf =1 :AnalysisEnd
    
    ratioTmp = ratioFurrow1Merged{lauf}
    ratioTmpFlank1 = ratioFurrow1{lauf}
    ratioTmpFlank2 = ratioFurrow2{lauf}
    
    ratioNorm{lauf} = (ratioTmp- ratioFurrowMean)/(ratioPoleMean- ratioFurrowMean)
    ratioNormFlank1{lauf} = (ratioTmpFlank1- ratioFurrowMean1)/(ratioPoleMean1- ratioFurrowMean1)
    ratioNormFlank2{lauf} = (ratioTmpFlank2- ratioFurrowMean2)/(ratioPoleMean2- ratioFurrowMean2)
    
end






%%%%%% The following section shifts the arrays of the linescan in a way
%%%%%% that every linescan has 300px in lenght and is centered. 
%%%%%% For the flanks
for lauf = 1: length(midOut)
    if midOut(lauf) < 150


        linescanSelectorTmp_I =  ratioNorm{lauf}

       % length(linescanSelectorTmp_I);
       % midOut(lauf);

        addZero = zeros(1,150 - midOut(lauf))';

        linescanSelectorCrop = cat(1,addZero,linescanSelectorTmp_I);
        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorCrop,addZero);

        linescanSelectorCropStoreGreen(:,lauf) = linescanSelectorCrop(1:300);
    elseif midOut(lauf) > 150

         linescanSelectorTmp_I = ratioNorm{lauf}


        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreen(:,lauf) =  linescanSelectorCrop(midOut(lauf)-149:300+(midOut(lauf)-150));


    elseif midOut(lauf) == 150
        
         linescanSelectorTmp_I =  ratioNorm{lauf}

        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1, linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreen(:,lauf) = linescanSelectorCrop(1:300);

    end
end




%%%%%% The following section shifts the arrays of the linescan in a way
%%%%%% that every linescan has 300px in lenght and is centered.
%%%%%% For the poles

for lauf = 1: length(midOutPole)
    if midOutPole(lauf) < 150


        linescanSelectorTmp_I =  ratioPoleMerged{lauf}

       % length(linescanSelectorTmp_I);
       % midOut(lauf);

        addZero = zeros(1,150 - midOutPole(lauf))';

        linescanSelectorCrop = cat(1,addZero,linescanSelectorTmp_I);
        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorCrop,addZero);

        linescanSelectorCropStoreGreenPole(:,lauf) = linescanSelectorCrop(1:300);
    elseif midOutPole(lauf) > 150

         linescanSelectorTmp_I = ratioPoleMerged{lauf}


        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreenPole(:,lauf) =  linescanSelectorCrop(midOutPole(lauf)-149:300+(midOutPole(lauf)-150));


    elseif midOutPole(lauf) == 150

        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1, linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreenPole(:,lauf) = linescanSelectorCrop(1:300);

    end
end


%%%%%% The following section shifts the arrays of the linescan in a way
%%%%%% that every linescan has 300px in lenght and is centered.

longDistanceVec = ((-149:150)*voxelX_mum)'

imshow(linescanSelectorCropStoreGreen,[])


colormap(jet)



normalizedKymo = normalizedImage(linescanSelectorCropStoreGreen);

imshow(normalizedKymo)


imwrite(normalizedKymo,[tifFilename,'Aniso_Kymo_Nrom.tif']);



 clear linescanSelectorTmp
    linescanSelectorTmp = zeros(300,43);
     linescanSelectorTmp  = double(linescanSelectorTmp)
midOut


linescanSelectorIn = redOut;

%%%%% This function generates the overlay
[greenOverlayStore]= doGenerateAnisotropieOverlay(greenNorm,  ratioNormFlank1,  ratioNormFlank2,flank1Store,flank2Store,AnalysisEnd);

  
         tiffwrite_mat(greenOverlayStore,[tifFilename,'Aniso_OverlayNorm.tif']);
          
         
         
         for lauf = 1:AnalysisEnd
             greenOverlayStoreTmp = greenOverlayStore(:,:,lauf).*(greenOverlayStore(:,:,lauf)>0);
             greenOverlayStoreNorm(:,:,lauf) =   ((greenOverlayStoreTmp - min(min(greenOverlayStoreTmp))  ./ (max(max(greenOverlayStoreTmp))- min(min(greenOverlayStoreTmp)))).*255);
       
         end
                 
         
         greenOverlayStoreNorm = normalizedImage3D(greenOverlayStoreNorm);
         
        tiffwrite_mat(greenOverlayStoreNorm,[tifFilename,'Aniso_OverlayNorm.tif']);
%        

       for lauf = 1:AnalysisEnd
        
       plot(1:length(linescanSelectorCropStoreGreen(:,lauf)),linescanSelectorCropStoreGreen(:,lauf))
       hold on
       
       end

       lauf = 1;
       
       %%%%%%% Calculate the the vector in the spline
I1_coords = flank1Store{lauf}';
I2_coords = flank2Store{lauf}';



spline = cscvn(I1_coords);
fnplt(spline,'or',2.5); hold on
t = 1:1:length(spline.breaks);
cv = fnval(spline, t);
cdv = fnval(fnder(spline), t);
quiver(cv(1,:),cv(2,:), cdv(1,:),cdv(2,:));


spline = cscvn(I2_coords);
fnplt(spline,'or',2.5); hold on
t = 1:1:length(spline.breaks);
cv = fnval(spline, t);
cdv = fnval(fnder(spline), t);
quiver(cv(1,:),cv(2,:), cdv(1,:),cdv(2,:),1);


%Generate quiver plot with scale = 1 and save as a lossless tiff

print('-dtiffn', 'quiver_cell9_frame5');

  for subrun= 1:length(cdv(1,:))
      
     
      
   v = sqrt((cdv(1,subrun))^2 + cdv(2,subrun)^2);

    alphaDeg(subrun) = acosd(cdv(1,subrun)/v);
   % IntenistyFlank1(subrun) = greenOverlayStoreNorm(round(cv(1,subrun)),round(cv(2,subrun)))
     
  end

  

%%%% Calculate the angle for the vector



  imshow(greenOverlayStoreNorm(:,:,lauf),[])
  hold on

  
  %%%% Correlate the tanges angle of the spline with each corresponding
  %%%% pixel
  flank1Tmp = flank1Store{lauf}';
  
  ratioNormFlank1Tmp = ratioNormFlank1{16};
   ratioNormFlank2Tmp = ratioNormFlank2{16}
  flankCrop1 = ratioNormFlank1Tmp(1:length(ratioNormFlank1Tmp(1,:)-5));
  flankCrop2 = ratioNormFlank1Tmp(1:length(ratioNormFlank2Tmp(1,:)-5));
  
%plot(alphaDeg1New,ratioNormFlank1Tmp(1:length(ratioNormFlank1Tmp(1,:)-5)),'x')

%plot(alphaDeg1New,ratioNormFlank2Tmp(1:length(ratioNormFlank1Tmp(1,:)-5)),'x')



%%%%%% Generate the Fit model

    I1_coords = flank1Store{16}';
    I2_coords = flank2Store{16}';
   
    
    %%%%%%%%% Esay angle detection method by simply calcualting the angle
    %%%%%%%%% between consecuteive points
    
    for lauf = 1:(length(I1_coords))-5
        
        uNew(lauf) = I1_coords(1,lauf+5) - I1_coords(1,lauf);
        vNew(lauf) = I1_coords(2,lauf+5) - I1_coords(2,lauf);
    
    end
    
    
    for subrun = 1:length(vNew)
    
            v1New = sqrt((vNew(1,subrun))^2 + uNew(1,subrun)^2);
            alphaDegModel(subrun) = acosd(uNew(1,subrun)/v1New)
                
    end
    
    %quiver(vNew(1,:),vNew(2,:), uNew(1,:),uNew(2,:),1);

[fitresult, gof] = createFitPolyFitAngleIntensity(alphaDegModel, ratioNormFlank1Tmp(1:end-2))

y180 = feval(fitresult,180)

y0 = feval(fitresult,0)

maxMean = (y180+y0)/2

%%%%% Correct Intensities

   expV_F1_Sammel = {}
   expV_F1_Sammel = {}
   f1_Diff = {}
   f2_Diff = {}
   alpha_1 = {}
   alpha_2  ={}

for lauf = 1:AnalysisEnd

    
  %%%% Correlate tanges angle with each pixel
  flank1Tmp = flank1Store{lauf}';
  flank2Tmp = flank2Store{lauf}';
  
  ratioNormFlank1Tmp = ratioNormFlank1{lauf};
  ratioNormFlank2Tmp = ratioNormFlank2{lauf};
  
  flankCrop1 = ratioNormFlank1Tmp(1:length(ratioNormFlank1Tmp(:,1)-5));
  flankCrop2 = ratioNormFlank2Tmp(1:length(ratioNormFlank2Tmp(:,1)-5));
  
    
    
    [expectedValues_F1,flank1DiffFromExpected_F1,alphaDeg_F1] =  doCalculateExpectedIntenisties(flank1Tmp,ratioNormFlank1Tmp,fitresult,flankCrop1);
    [expectedValues_F2,flank1DiffFromExpected_F2,alphaDeg_F2] =  doCalculateExpectedIntenisties(flank2Tmp,ratioNormFlank2Tmp,fitresult,flankCrop2);


   expV_F1_Sammel{lauf} = expectedValues_F1
   expV_F2_Sammel{lauf} = expectedValues_F2
   f1_Diff{lauf} = flank1DiffFromExpected_F1
   f2_Diff{lauf} = flank1DiffFromExpected_F2
   alpha_1{lauf} = alphaDeg_F1
   alpha_2{lauf}  = alphaDeg_F2
    
        
end

lauf = 9


[expOut midOut] = doAverageFlanksAngleCorrected(redMax1,redMax2,expV_F1_Sammel,expV_F2_Sammel)
[ratioMergeOut midOut] = doAverageFlanksAngleCorrected(redMax1,redMax2, ratioNormFlank1, ratioNormFlank2)


mkdir('predictedanisotropy')
curdir = pwd;

for lauf =1 :AnalysisEnd
    
  
     
    
    expMerge = expOut{lauf}
  ratioMerge =ratioMergeOut{lauf}
     
maxDist = round(length(expMerge)/2)
newDistanceVec = (0 - (maxDist-1):maxDist-1) *voxelX_mum


    lengthA = length(newDistanceVec);
    lengthB = length(expMerge);
    
    lengthVec = [lengthA lengthB];
    
    lengthEnd = min(lengthVec)
    
      
 %% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData(newDistanceVec(1:lengthEnd)',  expMerge(1:lengthEnd));

% Set up fittype and options.
ft = fittype( 'poly6' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf];
opts.Robust = 'LAR';
opts.Upper = [Inf Inf Inf Inf Inf Inf Inf];

% Fit model to data.
[fitOut, gofOut] = fit( xData, yData, ft, opts );

    
 h = figure(lauf)
 %plot(1:length(expectedValues_F1), diff_F1 ,'b');hold on
 plot(newDistanceVec,  expMerge(1:lengthEnd),'r');hold on
 plot(newDistanceVec, ratioMerge(1:lengthEnd),'b')
 plot(fitOut)
 
 
  
        
          set(gca,'FontSize',16,'FontName', 'Arial')
     % plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
       axis([-20 +20 0 1]) 
  title('Predictd and calculated ansiotropie at the cleavage furorw Exp. 2433','FontSize', 16,'FontName', 'Arial')
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial');
    ylabel('Normalized intensity  [A.U.]' ,'FontSize', 16,'FontName', 'Arial');
    legend(['Expected ratio at:' num2str(timeVec(lauf)), 's'] , ['Calculated ratio at:' num2str(timeVec(lauf)), 's'])
        pause(0.2)
 
        print(h,'-dpdf', [curdir '/' ,'predictedanisotropy' '/' , 'frame_', num2str(lauf) ,'_GaussFit_FWHM_anisotropy.pdf']);%tifCurvetifFilename);
    
end


figure(4)
plot(alphaTmp, correctedFlankTmp, 'x')
   
 
 figure(5)
 plot(alphaTmp,  correctedFlank1, 'x')


figure(3)
plot(alphaTmp, correctedFlank1,'x')




%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( alphaDeg, flankCrop );

       
    %%%%%%%%% Normalize data by subtracting 5% of flanking pixels first
    %%%%%%%%% selecting a 15µm window around the central pixel
    
   
    
    
    
    windowBoundaries = round(8/voxelX_mum);
    
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
    
   timeVec = timeVecTmp - (timeInterval * anaOnset)
   timeVec  =timeVecTmp;

%%%% Use the average of the total intensity for normalization
%%%% Has to be implemented
for lauf = 20:length(linescanSelectorIn)
 

   % lauf = 20
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

  % h=figure(1)
    %%%% Plot the data
  % plot(((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum),yGreenTmp,'c')
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
    
    for lauf = 20:AnalysisEnd
              h = figure(1) 
        plot((0-(newDistance-1):newDistance)*voxelX_mum,linescanSelector{lauf}')
       
     
      
        hold on
          set(gca,'FontSize',16,'FontName', 'Arial')
      plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
       axis([-15 +15 -1 1]) 
  title('Anisotropy SiR-actin green/red Exp. 2433','FontSize', 16,'FontName', 'Arial')
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial');
    ylabel('Ratio green/red  [A.U.]' ,'FontSize', 16,'FontName', 'Arial');
    legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(timeVec(lauf)) 's' ])
        pause(0.2)
       
       
      % print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5percAVG.pdf']);%tifCurvetifFilename);
   
      print(h,'-dpdf', [curdir '\' ,'redFit' '\' , 'frame_', num2str(lauf) ,'_GaussFit_FWHM_anisotropy.pdf']);%tifCurvetifFilename);
   close all
    end
    
    
    
         
    for lauf = 1:AnalysisEnd
              h = figure(1) 
        plot(((0-149):150)*voxelX_mum, linescanSelectorCropStoreGreen(:,lauf)')
       
     
      
        hold on
          set(gca,'FontSize',16,'FontName', 'Arial')
     % plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
       axis([-25 +25 0 1.5]) 
  title('anisotropy SiR-actin green/red Exp. 2433','FontSize', 16,'FontName', 'Arial')
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial');
    ylabel('Ratio green/red  [A.U.]' ,'FontSize', 16,'FontName', 'Arial');
    legend( [num2str(timeVec(lauf)), 's'] )
        pause(0.2)
       
       
      % print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5percAVG.pdf']);%tifCurvetifFilename);
   
      print(h,'-dpdf', [curdir '\' ,'redFit' '\' , 'frame_', num2str(lauf) ,'_FullLinescanLength_anisotropy.pdf']);%tifCurvetifFilename);
   close all
    end
    