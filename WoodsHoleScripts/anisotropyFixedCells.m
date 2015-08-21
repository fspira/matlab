%%%%%% This program will calculate the anisotropie at the equator of fixed
%%%%%% cell. The structure is to load the mid sections of the anisotropie
%%%%%% and the conventional fluorescence image, the user has to click on
%%%%%% the ingressing furrow and the anisotropie within a sliding window
%%%%%% along the cortex will be measured. The anisotropie at the poles will
%%%%%% be measured with a larger window than the one at the equator.

%% Header

clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;
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
%% Load images
curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell1_aniso.tif';
tifFilenameMid = 'cell1_zStack_Dapi-1.tif';

saveFileName = 'cell1_analysis';

%load('ratioAnisoParameters.mat')



%%%%% Load images

imgMidTmp = tiffread30(char(tifFilenameMid));

imgMidTmp = cat(3,imgMidTmp.data);
imgMid = imgMidTmp;
clear imgMidTmp

%%%%%% section required for multi stack images

voxelX_mumMid = 80; %size in nm

imgMidBlue = imgMid(:,:,1);
imgMidGreen = imgMid(:,:,2);


imgOrig = tiffread30(char(tifFilename))
imgOrigTmp = cat(3,imgOrig.data);

imgAniso = imgOrigTmp(:,:,1);
imgIntensity = imgOrigTmp(:,:,2);

voxelX_mum = 80; % pixel in nm

clear imagOrig

truncName = findstr(tifFilename,'.lsm');
emptyFlag = isempty(truncName);

if emptyFlag ==1

 truncName = findstr(tifFilename,'.tif');

end


folderName = tifFilename(1:truncName-1);
mkdir([curdir '/' folderName]);

cd([curdir '/' folderName]);

[m n p] = size(imgIntensity);

p=1

%% SEGMENT the images

t3Store = doSegmentImage(imgIntensity);

 
[m n p] = size(t3Store)
flank1Store = {};
flank2Store = {};
flank2Mid = [];
flank1Mid = [];

%%%%%% Detect Flanks

 
counter = 1;
for lauf = 1:p
    

    [flank1 flank2 flankMid1 flankMid2] = doSplitWatershedInTwo(t3Store(:,:,lauf), imgIntensity(:,:,lauf))

    flank1Store{lauf} = flank1
    flank2Store{lauf} = flank2
    flank1Mid(lauf) = flankMid1
    flank2Mid(lauf) = flankMid2
     counter = counter+1;
end
    
%%%%%% Burn ROI's into the image
%flank1Tmp = flank1Store{5};
%flank2Tmp = flank2Store{5};
%flank1Tmp = pole1Store{5};
%flank1Tmp = pole2Store{5};

%templateImage = S1(:,:,5);
%templateImage(:,:) = 0;
    
%for lauf = 1:length(flank1Tmp)
  
 %   templateImage(flank1Tmp(lauf,2),flank1Tmp(lauf,1)) = 255;
    
%end



%imwrite(templateImage,'splitSegment_cell9_Frame5_Pole2.tif');
%counter = 32
for lauf = counter:1
                
     flank1Store{lauf} = [50 50]; 
    flank2Store{lauf} = [50 50];
    flank1Mid(lauf) = 1
    flank2Mid(lauf) = 1
    
end


%%%%% Detect Poles

counter = 1;
for lauf = 1:1
    

    [pole1 pole2 poleMid1 poleMid2] = doSplitWatershedInTwoPoles(t3Store(:,:,lauf),imgIntensity(:,:,lauf))

    pole1Store{lauf} = pole1
    pole2Store{lauf} = pole2
    pole1Mid(lauf) = poleMid1
    pole2Mid(lauf) = poleMid2
    counter = counter+1;
end

 AnalysisEnd = 1;

%%%%% get new Values for the linsecan at the equator
for lauf =1:1
   
    
    MIJ.createImage(imgIntensity(:,:,lauf));
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
    
        
     MIJ.createImage(imgIntensity(:,:,lauf));
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








%%%%% This part orderes the segmented poles
for lauf = 1:1 %length(imgGauss)
     pole1Tmp = pole1Store{lauf};
    
    pole1FirstFrame = pole1Store{1};

    pole2FirstFrame = pole2Store{1};
    
    pole1FirstCenter = pole1Mid(1);
    pole2FirstCenter = pole2Mid(1);
    
    
   
    
    pole2Tmp = pole2Store{lauf};
    
    MIJ.createImage(imgIntensity(:,:,lauf));
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
    
     MIJ.createImage(imgIntensity(:,:,lauf));
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
    
               

                MIJ.createImage(imgIntensity(:,:,:,lauf));
                MIJ.run('8-bit');
              
               %MIJ.run('Stack to RGB');
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
                
                    imshow(imgIntensity(:,:,:,lauf))
                    hold on
                    plot(cNewMid1(lauf),rNewMid1(lauf),'xb')
                    plot(cNewMid2(lauf),rNewMid2(lauf),'xr')
                    pause(0.5)
                    
                end
                
         



%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%% Mark the poles
counter = 1
    
           for lauf = 1:1
     
      flank1Tmp = pole1Store{lauf};
      flank2Tmp = pole2Store{lauf};
    
    %%%%% Mark ingressing furrow
    
               

                MIJ.createImage(imgIntensity(:,:,:,lauf));
                MIJ.run('8-bit');
              
              % MIJ.run('Stack to RGB');
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
               
            
                   dist1 = pdist2([cChromatin1(lauf),rChromatin2(lauf)],[cChromatin2(lauf),rChromatin2(lauf)]);
                   
                   distChromatinChromatin(lauf) = dist1*voxelX_mum;
                   
                   
                 
              
  end
            
  
  %% Stage the cell measuring at the furrow and at the pole
  
  imshow(imgMid(:,:,3),[],'InitialMagnification',1000)

           
           [xx,yy] = ginput(6)
           
           
           close all
           
           pole1 = (pdist2([xx(1) yy(1)],[ xx(2) yy(2)])) *voxelX_mumMid
           pole2 = (pdist2([xx(3) yy(3)],[ xx(4) yy(4)])) *voxelX_mumMid
           
           meanPole = (pole1+pole2)/2
           
           cleavageFurrow = (pdist2([xx(5) yy(5)],[ xx(6) yy(6)])) *voxelX_mumMid
           
           poleFurrowRatio = cleavageFurrow/meanPole
           
           

%%   %%%%% Generate the line profil of the flanking regions
for lauf = 1:AnalysisEnd %length(redStack)
    
    
    
    %%%% Load XY coordinates of the two flanking regions into the tmp
    %%%% variable which is later used for analysis. This variable will be
    %%%% updated every loop
         
         flank1Tmp = flank1Store{lauf};
         flank2Tmp = flank2Store{lauf};
         
         %%%% Load the red image
         MIJ.createImage(imgIntensity(:,:,lauf));
         MIJ.run('setLine8');%%%% increase linescan width
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);%%%% Draw line into the image
         
         
         %%%%% Update InterpolatedLine
          clear logData
          MIJ.run('getSplineCoordinates')%%%% Get interpolated coordinates of the line
          logData = MIJ.getLog;%%%%% Reading MIJI log file to get coordinates
          
          MIJ.run('closeLogWindow')%%%% Close log file
         
          MIJ.run('getLinescanRed');%%%% Get intensities under the line
         yRed1{lauf} = MIJ.getColumn('y');%%%% Intensity values
         xRed1{lauf} = MIJ.getColumn('x');%%%% Running variable for each pixel can be discarded
        logDataTmp = char(logData);%%%%% convert values to characters - important to search for a string in the next step
        logDataTmp1 = strread(logDataTmp, '%s'); %%%% Sort log data using regular expression
        
         xLogTmp = logDataTmp1(2:4:end); %%%% Save X coordinates of the interpolated linescan
         yLogTmp = logDataTmp1(3:4:end);%%% Save Y coordinates of the interpolatd linescan
         
             
        xLog = zeros(length(xLogTmp),1,'double'); %%%% Pre allocate X and Y variables
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)%%%% Save XY coordinates of the interpolated linescan
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
         flank1StoreNew{lauf} = [xLog yLog]%%% Save coordinates into store array
       
         
         MIJ.createImage(imgIntensity(:,:,lauf));%%%% Repeat the same procedure for the second flanking region
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
           
            % flank1Tmp = flank1StoreNew{lauf};
         
             %%%%%% Get lineprofiles from the green channel
         MIJ.createImage(imgAniso(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreen1{lauf} = MIJ.getColumn('y');
         xGreen1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
            % flank2Tmp = flank2StoreNew{lauf};
         
         MIJ.createImage(imgAniso(:,:,lauf));
         MIJ.run('setLine8');
          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yGreen2{lauf} = MIJ.getColumn('y');
         xGreen2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
         
         
         
        
         
end

%%%%%% measure the poles
%%
    
for lauf = 1:AnalysisEnd %length(redStack)
    
    %%%%%% measure linescans for flanks
    
         flank1Tmp = pole1Store{lauf};
         flank2Tmp = pole2Store{lauf};
         
         
         MIJ.createImage(imgIntensity(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         
         
         %%%%% Update InterpolatedLine
          clear logData
          MIJ.run('getSplineCoordinates')
          logData = MIJ.getLog;
          
          MIJ.run('closeLogWindow')
         
          MIJ.run('getLinescanRed');
         yRed1Pole{lauf} = MIJ.getColumn('y');
         xRed1Pole{lauf} = MIJ.getColumn('x');
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
         
        pole1StoreNew{lauf} = [xLog yLog]
       
         
         MIJ.createImage(imgIntensity(:,:,lauf));
         MIJ.run('setLine8');
          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         
        
         
         %%%% Update interpolatedLine
              clear logData
              MIJ.run('getSplineCoordinates')
              logData = MIJ.getLog;

             MIJ.run('closeLogWindow')
          
            MIJ.run('getLinescanRed');
            yRed2Pole{lauf} = MIJ.getColumn('y');
            xRed2Pole{lauf} = MIJ.getColumn('x');
            
            
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

            pole2StoreNew{lauf} = [xLog yLog]
            
             MIJ.run('closeAllWindows');
         
         
            %%%%% Linescan at the poles
          
        pole1Tmp = pole1Store{lauf};
        pole2Tmp = pole2Store{lauf};
         
         
        MIJ.createImage(imgIntensity(:,:,lauf));
        MIJ.run('setLine8');
        MIJ.setRoi( [ pole1Tmp(:,1)'; pole1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
        MIJ.run('getLinescanRed');
        yRedPole1{lauf} = MIJ.getColumn('y');
        xRedPole1{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
        MIJ.createImage(imgIntensity(:,:,lauf));
        MIJ.run('setLine8');
         MIJ.setRoi( [ pole2Tmp(:,1)'; pole2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
        MIJ.run('getLinescanRed');
        yRedPole2{lauf} = MIJ.getColumn('y');
        xRedPole2{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
         
         
         
        MIJ.createImage(imgAniso(:,:,lauf));
        MIJ.run('setLine8');
         MIJ.setRoi( [ pole1Tmp(:,1)'; pole1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
        MIJ.run('getLinescanRed');
        yGreenPole1{lauf} = MIJ.getColumn('y');
        xGreenPole1{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
        MIJ.createImage(imgAniso(:,:,lauf));
        MIJ.run('setLine8');
         MIJ.setRoi( [ pole2Tmp(:,1)'; pole2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
        MIJ.run('getLinescanRed');
        yGreenPole2{lauf} = MIJ.getColumn('y');
        xGreenPole2{lauf} = MIJ.getColumn('x');
        MIJ.run('closeAllWindows');
         
        
         
end

%%

flank1Old = flank1Store;
flank2Old = flank2Store;

flank1Store = flank1StoreNew;
flank2Store = flank2StoreNew;


pole1Store = pole1StoreNew;
pole2Store = pole2StoreNew;

%%
%test = yRedMid1

%%%%%%% Identification of midsections of the ROI's. The following piece of
%%%%%%% code plots high intensity pixels at the center position into the
%%%%%%% linescan, which then can easily be detected by matlab.

%Miji;
%rMark1 = rMark1+1;
for lauf=1:AnalysisEnd
    %lauf =4
    redStackMid = imgIntensity(:,:,lauf);
    redStackMid(:,:) = 0;
    
    
         flank1Tmp = flank1Store{lauf};
         flank2Tmp = flank2Store{lauf};
         
    
    cMark1 = cNewMid1(lauf);
    rMark1 = rNewMid1(lauf);
    
    cMark2 =  cNewMid2(lauf);
    rMark2 =  rNewMid2(lauf);
     
    
         redStackMid(rMark1,cMark1) = 63000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid1{lauf} = MIJ.getColumn('y');
         xRedMid1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redMax1(lauf) = find(yRedMid1{lauf}==max(yRedMid1{lauf}))
         
         redStackMid(rMark2,cMark2) = 630000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid2{lauf} = MIJ.getColumn('y');
         xRedMid2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redMax2(lauf) = find(yRedMid2{lauf}==max(yRedMid2{lauf}))
         
         
         %%%%% Detection of the polar region
         
         
         pole1Tmp = pole1Store{lauf};
         pole2Tmp = pole2Store{lauf};
          
     
     cMark1 = cPoleMid1N(lauf);
     rMark1 = rPoleMid1N(lauf);
     
     cMark2 =  cPoleMid2N(lauf);
     rMark2 =  rPoleMid2N(lauf);
        redStackMid(:,:) = 0;
%     
          redStackMid(rMark1,cMark1) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ pole1Tmp(:,1)'; pole1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMidPole1{lauf} = MIJ.getColumn('y');
         xRedMidPole1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
%          
         redMaxPole1(lauf) = find(yRedMidPole1{lauf}==max(yRedMidPole1{lauf}))
%          
         redStackMid(rMark2,cMark2) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ pole2Tmp(:,1)'; pole2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMidPole2{lauf} = MIJ.getColumn('y');
         xRedMidPole2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
%          
        redMaxPole2(lauf) = find(yRedMidPole2{lauf}==max(yRedMidPole2{lauf}))
%          
%          
%          
%        

end

%%




close all

[anisoOut intensityOut midOut] = doAverageFlanks(redMax1,redMax2,yGreen1,yGreen2,yRed1,yRed2)


[anisoOutPole intensityOutPole midOutPole] = doAverageFlanks(redMaxPole1(1:1),redMaxPole2(1:1),yGreenPole1(1:1),yGreenPole2(1:1),yRedPole1(1:1),yRedPole2(1:1))

close all
plot(anisoOut{1})
hold on 
plot(intensityOut{1})
%% Plot the distance maps
distance = (1:length(intensityOut{1})) .* (voxelX_mum/1000)


   plot(earlyAna_perpendicular_lateral_OrientationCell,initialRecoilManualEarlyAnaPerpendicular_lateral,'bo','Markersize',12, 'Linewidth',2)
    
   
 plot(meta_perpendicular_lateral_OrientationCell,initialRecoilManualMetaPerpendicular_lateral,'bd','Markersize',12, 'Linewidth',2)

 %plot(proMeta_lateral_OrientationCell,initialRecoilManualProMeta_lateral,'bs','Markersize',12, 'Linewidth',2)

 
%% Plot intenisty at poles and cleavage furrow
fig = figure(2);
hold on
title({'Intensity and Anisotropy at the averaged cleavage furrow',...
             'at the confocal polscope in fixed cells'},'FontSize', 16)
              xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial')
[AX,H1,H2] = plotyy(distance,anisoOut{1}, distance,  intensityOut{1});

set(get(AX(2),'Ylabel'),'String','Intensity [A.U.]','FontSize', 16,'color','r','FontName', 'Arial') 
             %  set(get(AX(1),'Ylabel'),'String','Distance [µm]','FontSize', 16) 
 set(get(AX(1),'Ylabel'),'String','Anisotropy','FontSize', 16,'color','b','FontName', 'Arial') 
 %set(AX(1),'YLim',[0 2])
 % set(AX(2),'YLim',[0 22])
  set(AX,{'ycolor'},{'b';'r'})
set([H1(1)], 'color', 'b' );
set([H2(1)], 'color', 'r' );

set(fig, 'CurrentAxes', AX(2));
hold on;
plot(distance,intensityOut{1},'r','Markersize',12, 'Linewidth',2);

set(fig, 'CurrentAxes', AX(1));
hold on;
plot(distance,anisoOut{1},'Markersize',12, 'Linewidth',2);

         ;
                %ylabel ('Notch ratio [A.U.]','FontSize', 16);
                
        
            legend([H1 H2], ...
              'Anisotropy', ...
              'Intensity')%, ...
             % 'SiR actin notch ratio')%, 'MRLC', 'RhoA','alpha-Actinin' )%, ...
             % 'Chromosome-pole Distance')

                              
    print(fig,'-dpdf', ['ExpXXX_IntensityAnisotropyCleavageFurrow.pdf']);%tifCurvetifFilename);
    close all
    distancePole = (1:length(intensityOutPole{1})) .* (voxelX_mum/1000)
    
    fig = figure;
hold on
title({'Intensity and Anisotropy at the averaged poles',...
             'at the confocal polscope in fixed cells'},'FontSize', 16)
              xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial')
[AX,H1,H2] = plotyy(distancePole,anisoOutPole{1}, distancePole,  intensityOutPole{1});

set(get(AX(2),'Ylabel'),'String','Intensity [A.U.]','FontSize', 16,'color','r','FontName', 'Arial') 
             %  set(get(AX(1),'Ylabel'),'String','Distance [µm]','FontSize', 16) 
 set(get(AX(1),'Ylabel'),'String','Anisotropy','FontSize', 16,'color','b','FontName', 'Arial') 
 %set(AX(1),'YLim',[0 2])
 % set(AX(2),'YLim',[0 22])
  set(AX,{'ycolor'},{'b';'r'})
set([H1(1)], 'color', 'b' );
set([H2(1)], 'color', 'r' );

set(fig, 'CurrentAxes', AX(2));
hold on;
plot(distancePole,intensityOutPole{1},'r','Markersize',12, 'Linewidth',2);

set(fig, 'CurrentAxes', AX(1));
hold on;
plot(distancePole,anisoOutPole{1},'Markersize',12, 'Linewidth',2);

         ;
                %ylabel ('Notch ratio [A.U.]','FontSize', 16);
                
        
            legend([H1 H2], ...
              'Anisotropy', ...
              'Intensity')%, ...
             % 'SiR actin notch ratio')%, 'MRLC', 'RhoA','alpha-Actinin' )%, ...
             % 'Chromosome-pole Distance')

                              
    print(fig,'-dpdf', ['ExpXXX_IntensityAnisotropyPole.pdf']);%tifCurvetifFilename);
                
                
    
    
   %% Select sliding window along to measure anisotropy and Intensity 5µ large window along the merged cortex
   
   midOut
   slidingWindowLength = round(2000 / 80) % distance in nm
   SW = round(slidingWindowLength/2);
   
   intensityShape = intensityOut{1};
   anisoShape = anisoOut{1};
   poleIntensityShape = intensityOutPole{1};
   poleAnisoShape = anisoOutPole{1}
   
   intensityWindow = mean(intensityShape(midOut-(SW-1):midOut+SW));
   
   anisoWindow = mean(anisoShape(midOut-(SW-1):midOut+SW));
   
   poleIntensityWindow = mean(poleIntensityShape)
   
   poleAnisoWindow = mean(poleAnisoShape)
   
   ratioAnisoFurrowPole = anisoWindow/poleAnisoWindow
%%
   %%%% Write file to disk
    save([saveFileName,'.mat'])
   
            dummyVariable = (1:p)';
            dummyVariable(:) = 0;
           

            saveVariables = {};

            saveVariables{1} = poleFurrowRatio
            saveVariables{2} =  intensityWindow;
            
            saveVariables{3} =   poleIntensityWindow;
            saveVariables{4} =    anisoWindow;
            saveVariables{5} =   poleAnisoWindow;
            
          
            saveVariables{6} =  ratioAnisoFurrowPole;
           
            
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['stagingParameter, intensityCleavageFurrow, intensityPole, intensityFurrowPoleRatio',...
                'anisoCleavageFurrow, anisoPole, anisoFurrowPoleRatio']
            
            outid = fopen([saveFileName,'Analysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([saveFileName,'Analysis.csv'],csvData,'roffset',1,'-append')
            
            
             cd([curdir]) 
  
            
