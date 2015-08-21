%%%%% anisotropy Linescan
%%
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
addpath('/Users/spira/Desktop/programme/WoodsHoleScripts')

javaaddpath '/Applications/Fiji.app/jars/ij-1.49h.jar'
addpath('/Applications/Fiji.app/scripts')
addpath('/Applications/Fiji.app/plugins/')
addpath('/Users/spira/Desktop/programme/WoodsHoleScripts/matlabScriptsShalin')
%addpath('/Users/spira/Desktop/2015_WoodsHole/matlabScriptsShalin')

%addpath('/Users/spira/Desktop/programme/miji/mij.jar')
%addpath('/Users/spira/Desktop/programme/miji/ij.jar')

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/mij.jar'

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/ij.jar'


 MIJ.start('/Applications/Fiji.app')


curdir = pwd;

origDir = '/Users/spira/Desktop/5FramePoleScope/20150503_SirActin_5Frame_dish2/cell1'

fileList= dir([curdir '/' '*I0_*'])

isopath='/Users/spira/Desktop/5FramePoleScope/calib256';
%normFactors=computeConfocalNormalizations(isopath,633,1.4,110);


clear avgFile




for lauf = 1 : length(fileList)
 
   tmpFile = tiffread30(fileList(lauf).name);
   tmpFileFluor = tmpFile.data{1};
    
    avgFile(:,:,lauf) =   tmpFileFluor;
   
end

tiffwrite_mat(avgFile, [origDir '/' 'analysis' '/' 'I0Channel.tif'])

cd(origDir)

fileListOrig = dir([origDir '/'  '*I*'])



  
[m n p] =size(avgFile)

AnalysisEnd =p
firstFrame =1
lastFrame = 22;
anaOnset =1;


saveFileName = 'cell1Analysis';

%%% load images

lengthDirectory = length(fileListOrig)






fileListOrig0 = dir([origDir '/'  '*I0_*'])
fileListOrig45 = dir([origDir '/'  '*I45_*'])
fileListOrig90 =  dir([origDir '/'  '*I90_*'])
fileListOrig135 = dir([origDir '/'  '*I135_*'])
fileListOrigB = dir([origDir '/'  '*I0Bleach_*'])
fileListOrigH2B = dir([origDir '/' '*H2B_*'])



          
  
p = length(fileListOrig0)

AnalysisEnd =p
%firstFrame = 22;
%lastFrame = 42;
%anaOnset =5;


%%% load images

lengthDirectory = length(fileListOrig0)

            clear listDirectory sortOut
            clear sortedList sortArray sortIdx

            for subrun = 1:length(fileListOrig0)-1
                laufName = fileListOrig0(subrun+1).name
                listDirectory{subrun} =  laufName
                suffixStrSort = findstr( laufName,'.lsm');
                prefixStrSort = findstr( laufName,'_');


                sortArray(subrun) = str2num(laufName(prefixStrSort+1:suffixStrSort-1));
            end

            [sortArray0 listDirectory0] = doSortDirectory(fileListOrig0)
            [sortOut0 sortIdx0] = sort(sortArray0)
           
            [sortArray45 listDirectory45] = doSortDirectory(fileListOrig45)
            [sortOut45 sortIdx45] = sort(sortArray45)
           
            [sortArray90 listDirectory90] = doSortDirectory(fileListOrig90)
            [sortOut90 sortIdx90] = sort(sortArray90)
           
            [sortArray135 listDirectory135] = doSortDirectory(fileListOrig135)
            [sortOut135 sortIdx135] = sort(sortArray135)
           
            [sortArrayB listDirectoryB] = doSortDirectory(fileListOrigB)
            [sortOutB sortIdxB] = sort(sortArrayB)
           
            [sortArrayH2B listDirectoryH2B] = doSortDirectory(fileListOrigH2B)
            [sortOutH2B sortIdxH2B] = sort(sortArrayH2B)
           

            sortedList0{1} = fileListOrig0(1).name
            sortedList45{1} = fileListOrig45(1).name
            sortedList90{1} = fileListOrig90(1).name
            sortedList135{1} = fileListOrig135(1).name
            sortedListB{1} = fileListOrigB(1).name
            sortedListH2B{1} = fileListOrigH2B(1).name

            for subrun = 1:length(sortIdx0)
                sortedList0{subrun+1} = listDirectory0{sortIdx0(subrun)};
                sortedList45{subrun+1} = listDirectory45{sortIdx45(subrun)};
                sortedList90{subrun+1} = listDirectory90{sortIdx90(subrun)};
                sortedList135{subrun+1} = listDirectory135{sortIdx135(subrun)};
                sortedListB{subrun+1} = listDirectoryB{sortIdxB(subrun)};
                sortedListH2B{subrun+1} = listDirectoryH2B{sortIdxH2B(subrun)};
            end

            
            Frame = length(fileListOrig0);
            p = length(fileListOrig0);
            
            
            
  
for lauf = 1:p
    I0Name(:,:,:,lauf) = tiffread30(char(sortedList0((lauf))));
    I45Name(:,:,:,lauf) = tiffread30(char(sortedList45((lauf))));
    I90Name(:,:,:,lauf) = tiffread30(char(sortedList90((lauf))));
    I135Name(:,:,:,lauf) = tiffread30(char(sortedList135((lauf))));
    IBName(:,:,:,lauf) = tiffread30(char(sortedListB((lauf))));
    IH2BName(:,:,:,lauf) = tiffread30(char(sortedListH2B((lauf))));

end

for lauf = 1:p
 
    I0Tmp = I0Name(:,:,:,lauf);
    I0Tmp = cat(3,I0Tmp.data);
    I0File(:,:,lauf) =  I0Tmp{1};
    
    I45Tmp = I45Name(:,:,:,lauf);
    I45Tmp = cat(3,I45Tmp.data);
    I45File(:,:,lauf) = I45Tmp{1};
    
    I90Tmp = I90Name(:,:,:,lauf);
    I90Tmp = cat(3,I90Tmp.data);
    I90File(:,:,lauf) = I90Tmp{1};
    
    I135Tmp = I135Name(:,:,lauf);
    I135Tmp = cat(3,I135Tmp.data);
    I135File(:,:,lauf) = I135Tmp{1};
    
     IBTmp = IBName(:,:,lauf);
    IBTmp = cat(3,IBTmp.data);
    IBFile(:,:,lauf) = IBTmp{1};
    
    
    IH2BTmp = IH2BName(:,:,lauf);
    IH2BTmp = cat(3,IH2BTmp.data);
    IH2BFile(:,:,lauf) = IH2BTmp{1};
 

end


I0Segment = doSegmentImage(I0File);


I45Segment = doSegmentImage(I45File);

I90Segment = doSegmentImage(I90File);
I135Segment = doSegmentImage(I135File);
IBSegment = doSegmentImage(IBFile);


for lauf = 1:p

    averageImage(:,:,lauf)= (I0File(:,:,lauf) + I45File(:,:,lauf)+ I90File(:,:,lauf)...
         + I135File(:,:,lauf)+ IBFile(:,:,lauf))/5;

end







for lauf = 1:p
   
    
     greenStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(averageImage(:,:,lauf)),green);
     redStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(IH2BFile(:,:,lauf)),red);
      blueStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(IH2BFile(:,:,lauf)),blue);
     imgMerge(:,:,:,lauf) = greenStackRGB(:,:,:,lauf) + redStackRGB(:,:,:,lauf)+blueStackRGB(:,:,:,lauf);
end


tiffwrite_mat(averageImage, [origDir '/' 'analysis' '/' 'averageImage.tif'])

tiffwrite_mat(I0File, [origDir '/' 'analysis' '/' 'I0File.tif'])
tiffwrite_mat(I45File, [origDir '/' 'analysis' '/' 'I45File.tif'])
tiffwrite_mat(I90File, [origDir '/' 'analysis' '/' 'I90File.tif'])
tiffwrite_mat(I135File, [origDir '/' 'analysis' '/' 'I135File.tif'])
tiffwrite_mat(IBFile, [origDir '/' 'analysis' '/' 'IBleachFile.tif'])
tiffwrite_mat(IH2BFile, [origDir '/' 'analysis' '/' 'H2BFile.tif'])


tiffwrite_RBG(imgMerge, [origDir '/' 'analysis' '/' 'imgMerge.tif'])


%imshow(averageImage(:,:,1),[]);



voxelX = getfield(I0Name(:,:,:,lauf),'lsm','VoxelSizeX');
voxelX_mum = voxelX*1000000;

timeInterval = 30;

anaTime = anaOnset*timeInterval;
timeMax = (p*timeInterval+p) - anaTime;
anaTime = 0;
 
 timeVec = 1:timeInterval:p*timeInterval;
            
timeVec = timeVec-timeVec(anaOnset);


%voxelX_mum
%voxelX_H2B

%%%% Select first channel from file


%%%% Crop the file

    clear I0Tmp
    I0Tmp = I0File(:,:,firstFrame:lastFrame);
    clear I0File
    I0File = I0Tmp;
    clear I0Tmp
    
    clear I45Tmp
    I45Tmp = I45File(:,:,firstFrame:lastFrame);
    clear I45File
    I45File = I45Tmp;
    clear I45Tmp
    
    clear I90Tmp
    I90Tmp = I90File(:,:,firstFrame:lastFrame);
    clear I90File
    I90File = I90Tmp;
    clear I90Tmp
    
    clear I135Tmp
    I135Tmp = I135File(:,:,firstFrame:lastFrame);
    clear I135File
    I135File = I135Tmp;
    clear I135Tmp
    
 

%%

%%
%% Background correction I135File

[I0File I45File I90File I135File I135Fileb] =  bkgCorrection5FrameConfocalPolscope(I0File,  I45File,I90File, I135File, IBFile);
%%%% average files

%%%%%%% sample file
greenStack = double(averageImage);


%greenStack = greenImgStore;
%redStack = redImgStore;

greenImgStore = greenStack;

greenImg = greenStack;
[m n p] = size(greenImg(:,:,:));
shiftDistance = 512-400 %m;

greenImg(1,512,:) = 0;greenImg(512,1,:) = 0;

imgOrigGreenTranslate =[];



imgOrigGreenTranslate = zeros(512,512,1:AnalysisEnd);


for lauf = 1:AnalysisEnd

    imgOrigGreenTranslate(:,:,lauf) = imtranslate(greenImg(:,:,lauf),[shiftDistance, shiftDistance]);
  

end


segment = doSegmentImage( imgOrigGreenTranslate);


%Miji;

for lauf = 1:p
 
    %%%%% This piece of code determines the axis of the cell and suggests
    %%%%% position of pole-pole and furrow-furrow centers
   % lauf = 15
  L2 = segment(:,:,lauf);
  cc = regionprops(L2,'all');
  B = bwboundaries(L2,8); 
  boundary = B{1};

  bStore{lauf} = B{1};
  
k=1
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
    
xbar = cc(k).Centroid(1);
    ybar = cc(k).Centroid(2);

    a = cc(k).MajorAxisLength/2;
    b = cc(k).MinorAxisLength/2;

    theta = pi*(cc(k).Orientation/180);
    
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;

    imshow(L2,[])
    hold on
    plot(x,y,'r','LineWidth',2);


   % tan((cc.Orientation+90)/360*2*pi)*250
   
   %%%%% left
   
   
   x = cc.Centroid(1,1)-1
    y = cc.Centroid(1,2)-(tan((cc.Orientation+90)/360*2*pi)*-1)
    plot(x,y,'dw')
     
    setLength = 1
    
      while pdist2([cc.Centroid(1,1),cc.Centroid(1,2)],[x,y]) < 180
      
          setLength = setLength +1
           x = cc.Centroid(1,1)-setLength
           y = cc.Centroid(1,2)-(tan((cc.Orientation+90)/360*2*pi)*-setLength)
  
      end
          
    
    [xCoord,yCoord,regionIntersect,selectedRoi] = detectAxis_v1(L2, cc,x,y,boundary,400)
    cNewMid1(lauf) = xCoord;
    rNewMid1(lauf) = yCoord;
    
      flank1 = selectedRoi;
      flankMid1 = regionIntersect;
      
    
   %%%%% right
    setLength =1;
    x = cc.Centroid(1,1)+1
    y = cc.Centroid(1,2)-(tan((cc.Orientation-90)/360*2*pi)*1)
    plot(x,y,'xg')
    
    
      
    setLength = 1
    
      while pdist2([cc.Centroid(1,1),cc.Centroid(1,2)],[x,y]) < 180
      
          setLength = setLength +1
           x = cc.Centroid(1,1)+setLength
           y = cc.Centroid(1,2)-(tan((cc.Orientation-90)/360*2*pi)*+setLength)
  
      end
    
    [xCoord, yCoord, regionIntersect,selectedRoi] = detectAxis_v1(L2, cc,x,y,boundary,400)
    cNewMid2(lauf) = xCoord;
    rNewMid2(lauf) = yCoord;
    
    
    
    flank2 = selectedRoi;
    flankMid2 = regionIntersect;
    
    %%%% lower
    
  
    
            setLength = 0.5;
    
          x = cc.Centroid(1,1)+setLength
          y = cc.Centroid(1,2)-(tan((cc.Orientation)/360*2*pi)*setLength)
  
      
        while pdist2([cc.Centroid(1,1),cc.Centroid(1,2)],[x,y]) <180
      
              setLength = setLength +0.5
            x = cc.Centroid(1,1)+setLength
            y = cc.Centroid(1,2)-(tan((cc.Orientation)/360*2*pi)*setLength)
  
        end
        
    
              %%% upper
        
                   [xCoord,yCoord,regionIntersect,selectedRoi] = detectAxis_v1(L2, cc,x,y,boundary,400)
         
                   cNewLongAxis1(lauf) =xCoord;
                   rNewLongAxis1(lauf) =yCoord;
                 
                      pole1 = selectedRoi;
                      poleMid1 = regionIntersect;
            
       
                %%%%% lower
                
                    setLength = 0.5;
    
          x = cc.Centroid(1,1)-setLength
          y = cc.Centroid(1,2)-(tan((cc.Orientation)/360*2*pi)*-setLength)
  
      
        while pdist2([cc.Centroid(1,1),cc.Centroid(1,2)],[x,y]) <180
      
              setLength = setLength +0.5
              x = cc.Centroid(1,1)-setLength
              y = cc.Centroid(1,2)-(tan((cc.Orientation)/360*2*pi)*-setLength)
         end
          
                  [xCoord,yCoord,regionIntersect,selectedRoi] = detectAxis_v1(L2, cc,x,y,boundary,400)
                  
                   cNewLongAxis2(lauf) =xCoord;
                   rNewLongAxis2(lauf) =yCoord;
                
                       pole2 = selectedRoi;
                       poleMid2 = regionIntersect;
        
        
             
        
        
   
   
     
      
      
      flank1Store{lauf} =flank1-shiftDistance;
      flank2Store{lauf} = flank2-shiftDistance;
      flankMid1Store{lauf} = flankMid1-shiftDistance;
      flankMid2Store{lauf}  = flankMid2-shiftDistance;
      
      
      pole1Store{lauf} = pole1-shiftDistance;
      pole2Store{lauf} = pole2-shiftDistance;
      
      poleMid1Store{lauf} = poleMid1-shiftDistance;
      poleMid2Store{lauf} = poleMid2-shiftDistance;
      
     
    
end



%%%%% Shift segmentation coordinates


     

fluoStart = 22
%clear cSort1 rSort1 cSort2 rSort2
  for lauf = 1:fluoStart-1  %%% if auto detection block is used, this goes to  fluoStart-1
    flank1StoreTmp = flank1Store{lauf};
    flank2StoreTmp = flank2Store{lauf};
    cSort1{lauf} =  flank1StoreTmp(:,2)
    
    rSort1{lauf}    = flank1StoreTmp(:,1)
    cSort2{lauf}    = flank2StoreTmp(:,2);
    rSort2{lauf}    = flank2StoreTmp(:,1);
 end
 
 for lauf = fluoStart:p  %%% if auto detection block is used, this goes to  fluoStart-1
    cSort1{lauf} = 50;
    rSort1{lauf} = 50;
    cSort2{lauf} = 50;
    rSort2{lauf} = 50;
 end
% clear cNonPol1 rNonPol2 cNonPol2 rNonPol2
  for lauf = 1:fluoStart-1  %%% if auto detection block is used, this goes to  fluoStart-1
   
   pole1StoreTmp = pole1Store{lauf};
   pole2StoreTmp = pole2Store{lauf};
      
   cNonPol1{lauf} = pole1StoreTmp(:,2);
   rNonPol1{lauf} = pole1StoreTmp(:,1);
   cNonPol2{lauf} = pole2StoreTmp(:,2);
   rNonPol2{lauf} = pole2StoreTmp(:,1);
   
 end
 
 for lauf = fluoStart:p  %%% if auto detection block is used, this goes to  fluoStart-1
    cNonPol1{lauf} = 50;
   rNonPol1{lauf}    = 50;
     cNonPol2{lauf}    = 50;
   rNonPol2{lauf}    = 50;
 end
 

 for lauf = fluoStart:p  %%% if auto detection block is used, this goes to  fluoStart-1
   cNewMid1(lauf) = 50;
   rNewMid1(lauf) =50;
   cNewLongAxis1(lauf) = 50;
   rNewLongAxis1(lauf) = 50;
   
   cNewMid2(lauf) = 50;
   rNewMid2(lauf) =50;
   cNewLongAxis2(lauf) = 50;
   rNewLongAxis2(lauf) = 50;
 end
 %Miji;
 
 [m n p] = size(greenStack);
 
 
 %%%% Draw linescan along the cleavage furrow
 
clear  flank1TmpUpdate flank2TmpUpdate cNew1Ref rNew1Ref...
     cNew2Ref rNew2Ref flank1StoreNew flank2StoreNew...
     flank1Store flank2Store

 
 
for lauf =1:p %length(imgGauss)
    
    
   clear flank1TmpUpdate flank2TmpUpdate
    
     MIJ.createImage(imgMerge(:,:,:,lauf));
      MIJ.run('8-bit');
    MIJ.run('Stack to RGB');
    MIJ.setRoi( [cSort1{lauf}';rSort1{lauf}'], ij.gui.Roi.POLYLINE);
   
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
    
     cNew1Ref{lauf} =  xLog;
     rNew1Ref{lauf} =  yLog;
    
    
    flank1Store{lauf} = [flank1TmpUpdate(:,1) flank1TmpUpdate(:,2)]
    %%%%%%%
    
        
     MIJ.createImage(imgMerge(:,:,:,lauf));
       MIJ.run('8-bit');
    MIJ.run('Stack to RGB');
     MIJ.setRoi( [cSort2{lauf}';rSort2{lauf}'], ij.gui.Roi.POLYLINE);
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
      cNew2Ref{lauf} =  xLog;
      rNew2Ref{lauf} =  yLog;
    
 %  MIJ.run('saveMacro');
 %   sROI= ReadImageJROI('Roi.zip');
 %   tmp = sROI{:,:}.mnCoordinates;
  %  flank2TmpUpdate(:,1) = tmp(:,1);
  %  flank2TmpUpdate(:,2) = tmp(:,2);
     % MIJ.closeAllWindows
    flank2Store{lauf} = [flank2TmpUpdate(:,1) flank2TmpUpdate(:,2)]
    
    cNew2Ref{lauf} =  xLog;
    rNew2Ref{lauf} =  yLog;
    
    MIJ.run('closeAllWindows');
  
    lauf
    
    
    
    

end


for lauf = 1:p %length(imgGauss)
    
    
    clear flank1TmpUpdate flank2TmpUpdate
    
    
    
    
    
    MIJ.createImage(imgMerge(:,:,:,lauf));
    MIJ.run('8-bit');
    MIJ.run('Stack to RGB');
    MIJ.setRoi( [ cNonPol1{lauf}'; rNonPol1{lauf}'], ij.gui.Roi.POLYLINE);
    
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
    
    
    pole1StoreNew{lauf} = [xLog yLog]
    
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    flank1TmpUpdate(:,1) = xLog;
    flank1TmpUpdate(:,2) = yLog;
    
    cNonPol1Store{lauf} =  xLog;
    rNonPol1Store{lauf} =  yLog;
    
    
    pole1Store{lauf} = [flank1TmpUpdate(:,1) flank1TmpUpdate(:,2)]
    %%%%%%%
    
    
    MIJ.createImage(imgMerge(:,:,:,lauf));
    MIJ.run('8-bit');
    MIJ.run('Stack to RGB');
    MIJ.setRoi( [ cNonPol2{lauf}'; rNonPol2{lauf}'], ij.gui.Roi.POLYLINE);
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
    
    pole2StoreNew{lauf} = [xLog yLog]
    
    
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    flank2TmpUpdate(:,1) = xLog;
    flank2TmpUpdate(:,2) = yLog;
    
    %  MIJ.run('saveMacro');
    %   sROI= ReadImageJROI('Roi.zip');
    %   tmp = sROI{:,:}.mnCoordinates;
    %  flank2TmpUpdate(:,1) = tmp(:,1);
    %  flank2TmpUpdate(:,2) = tmp(:,2);
    % MIJ.closeAllWindows
    pole2Store{lauf} = [flank2TmpUpdate(:,1) flank2TmpUpdate(:,2)]
    
    cNonPol2Store{lauf} =  xLog;
    rNonPol2Store{lauf} =  yLog;
    
    
    
    
    
    
    MIJ.run('closeAllWindows');
    lauf
    
end



 
%%%%%% Detect pixels which are most close together
clear dist1;
clear minIndex;
clear  distMinStore


            
            
           % cd(curdir);
         %    Miji;
         % lauf =2;
         for lauf = 1:p%length(imgGauss)
             
             
             %%%%% Mark ingressing furrow
             
             
             
             MIJ.createImage(imgMerge(:,:,:,lauf));
             MIJ.run('8-bit');
             
             MIJ.run('Stack to RGB');
             % MIJ.selectWindow('Import from Matlab');
             MIJ.setRoi( [cNewMid1(lauf);rNewMid1(lauf)], ij.gui.Roi.POINT);
             k = waitforbuttonpress
             coords =    MIJ.getRoi(1);
             cNewMid1N(lauf) = coords(1);
             rNewMid1N(lauf) = coords(2);
             
             
             
             MIJ.setRoi( [cNewMid2(lauf);rNewMid2(lauf)], ij.gui.Roi.POINT);
             k = waitforbuttonpress
             coords =    MIJ.getRoi(1);
             cNewMid2N(lauf) = coords(1);
             rNewMid2N(lauf) = coords(2);
             
             
             MIJ.setRoi( [cNewLongAxis1(lauf);rNewLongAxis1(lauf)], ij.gui.Roi.POINT);
             
             k = waitforbuttonpress
             
             coords =    MIJ.getRoi(1);
             
             cNewLongAxis1(lauf) = coords(1);
             rNewLongAxis1(lauf) = coords(2);
             
             
             MIJ.setRoi( [cNewLongAxis2(lauf);rNewLongAxis2(lauf)], ij.gui.Roi.POINT);
             
             k = waitforbuttonpress
             coords =    MIJ.getRoi(1);
             
             cNewLongAxis2(lauf) = coords(1);
             rNewLongAxis2(lauf) = coords(2);
             
             
             
             MIJ.run('closeAllWindows');
             
             
             
             
             
             
         end
         
         
              cNewMid1 =  cNewMid1N;
                rNewMid1 =  rNewMid1N;
                
                
                cNewMid2 =  cNewMid2N;
                rNewMid2 =  rNewMid2N;


                
                
for lauf = 1:AnalysisEnd %length(redStack)
        
    
    
    
     
    %%%% Load XY coordinates of the two flanking regions into the tmp
    %%%% variable which is later used for analysis. This variable will be
    %%%% updated every loop
         
         flank1Tmp = flank1StoreNew{lauf};
         flank2Tmp = flank2StoreNew{lauf};
         
         %%%% Load the red image
         MIJ.createImage(I0File(:,:,lauf));
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
       
         
         MIJ.createImage(I0File(:,:,lauf));%%%% Repeat the same procedure for the second flanking region
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
         MIJ.createImage(I0File(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI0_1{lauf} = MIJ.getColumn('y');
         xI0_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
              MIJ.createImage(I45File(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,2)'; flank1Tmp(:,1)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI45_1{lauf} = MIJ.getColumn('y');
         xI45_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
             MIJ.createImage(I90File(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI90_1{lauf} = MIJ.getColumn('y');
         xI90_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
              MIJ.createImage(I45File(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,2)'; flank1Tmp(:,1)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI135_1{lauf} = MIJ.getColumn('y');
         xI135_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
          
              MIJ.createImage(I45File(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,2)'; flank1Tmp(:,1)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yIB_1{lauf} = MIJ.getColumn('y');
         xIB_1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
        
         
         
         %%%%%% Flank2
         
          MIJ.createImage(I0File(:,:,lauf));
         MIJ.run('setLine8');
          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI0_2{lauf} = MIJ.getColumn('y');
         xI0_2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
              MIJ.createImage(I45File(:,:,lauf));
         MIJ.run('setLine8');
        MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI45_2{lauf} = MIJ.getColumn('y');
         xI45_2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
             MIJ.createImage(I90File(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI90_2{lauf} = MIJ.getColumn('y');
         xI90_2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
              MIJ.createImage(I45File(:,:,lauf));
         MIJ.run('setLine8');
        MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yI135_2{lauf} = MIJ.getColumn('y');
         xI135_2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
          
              MIJ.createImage(I45File(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yIB_2{lauf} = MIJ.getColumn('y');
         xIB_2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
         
%          
%          
%          %%%%%%%%%%%
%     
%     %%%%% Repeat the same anlaysis as for the flanks for the pole
%     
%     
%      %%%% Load XY coordinates of the two flanking regions into the tmp
%     %%%% variable which is later used for analysis. This variable will be
%     %%%% updated every loop
%          
%          flank1Tmp = pole1StoreNew{lauf};
%          flank2Tmp = pole2StoreNew{lauf};
%          
%          %%%% Load the red image
%          MIJ.createImage(redStack(:,:,lauf));
%          MIJ.run('setLine8');%%%% increase linescan width
%          MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);%%%% Draw line into the image
%          
%          
%          %%%%% Update InterpolatedLine
%           clear logData
%           MIJ.run('getSplineCoordinates')%%%% Get interpolated coordinates of the line
%           logData = MIJ.getLog;%%%%% Reading MIJI log file to get coordinates
%           
%           MIJ.run('closeLogWindow')%%%% Close log file
%          
%           MIJ.run('getLinescanRed');%%%% Get intensities under the line
%          yRed1{lauf} = MIJ.getColumn('y');%%%% Intensity values
%          xRed1{lauf} = MIJ.getColumn('x');%%%% Running variable for each pixel can be discarded
%         logDataTmp = char(logData);%%%%% convert values to characters - important to search for a string in the next step
%         logDataTmp1 = strread(logDataTmp, '%s'); %%%% Sort log data using regular expression
%         
%          xLogTmp = logDataTmp1(2:4:end); %%%% Save X coordinates of the interpolated linescan
%          yLogTmp = logDataTmp1(3:4:end);%%% Save Y coordinates of the interpolatd linescan
%          
%              
%         xLog = zeros(length(xLogTmp),1,'double'); %%%% Pre allocate X and Y variables
%         yLog = zeros(length(yLogTmp),1,'double');
%          
%          for subrun = 1 : length(xLogTmp)%%%% Save XY coordinates of the interpolated linescan
%              
%              xLog(subrun) = round(str2num(xLogTmp{subrun}));
%              yLog(subrun) = round(str2num(yLogTmp{subrun}));
%              
%          end
%          
%          pole1StoreNew{lauf} = [xLog yLog]%%% Save coordinates into store array
%        
%          
%          MIJ.createImage(redStack(:,:,lauf));%%%% Repeat the same procedure for the second flanking region
%          MIJ.run('setLine8');
%           MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%          
%         
%          
%          %%%% Update interpolatedLine
%               clear logData
%               MIJ.run('getSplineCoordinates')
%               logData = MIJ.getLog;
% 
%              MIJ.run('closeLogWindow')
%           
%             MIJ.run('getLinescanRed');
%             yRed2{lauf} = MIJ.getColumn('y');
%             xRed2{lauf} = MIJ.getColumn('x');
%             
%             
%             logDataTmp = char(logData);
%             logDataTmp1 = strread(logDataTmp, '%s');
% 
%              xLogTmp = logDataTmp1(2:4:end);
%              yLogTmp = logDataTmp1(3:4:end);
% 
% 
%             xLog = zeros(length(xLogTmp),1,'double');
%             yLog = zeros(length(yLogTmp),1,'double');
% 
%              for subrun = 1 : length(xLogTmp)
% 
%                  xLog(subrun) = round(str2num(xLogTmp{subrun}));
%                  yLog(subrun) = round(str2num(yLogTmp{subrun}));
% 
%              end
% 
%              pole2StoreNew{lauf} = [xLog yLog]
%             
%              MIJ.run('closeAllWindows');
%            
%             % flank1Tmp = flank1StoreNew{lauf};
%          
%              %%%%%% Get lineprofiles from the green channel
%          MIJ.createImage(greenStack(:,:,lauf));
%          MIJ.run('setLine8');
%          MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%          MIJ.run('getLinescanRed');
%          yGreenNonPol1{lauf} = MIJ.getColumn('y');
%          xValuesNonPol1{lauf} = MIJ.getColumn('x');
%          MIJ.run('closeAllWindows');
%          
%          
%          
%          
%          MIJ.createImage(redStack(:,:,lauf));
%          MIJ.run('setLine8');
%          MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%          MIJ.run('getLinescanRed');
%          yRedNonPol1{lauf} = MIJ.getColumn('y');
%          xRedNonPol1{lauf} = MIJ.getColumn('x');
%          MIJ.run('closeAllWindows');
%          
%          
%          
%             % flank2Tmp = flank2StoreNew{lauf};
%          
%          MIJ.createImage(greenStack(:,:,lauf));
%          MIJ.run('setLine8');
%           MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%          MIJ.run('getLinescanRed');
%          yGreenNonPol2{lauf} = MIJ.getColumn('y');
%          xValuesNonPol2{lauf} = MIJ.getColumn('x');
%          MIJ.run('closeAllWindows');
%          
%             MIJ.createImage(redStack(:,:,lauf));
%          MIJ.run('setLine8');
%          MIJ.setRoi( [ flank2Tmp(:,1)'; flank2Tmp(:,2)'], ij.gui.Roi.POLYLINE);
%          MIJ.run('getLinescanRed');
%          yRedNonPol2{lauf} = MIJ.getColumn('y');
%          xRedNonPol2{lauf} = MIJ.getColumn('x');
%          MIJ.run('closeAllWindows');
%          
%     
%     
    
    
    
    
         
         
         
end


figure(1)
hold on
for lauf = 1:33
plot(yI0_1{lauf})
end



for lauf=1:AnalysisEnd
    %lauf =4
       
         flank1Tmp = flank1StoreNew{lauf};
         flank2Tmp = flank2StoreNew{lauf};
        % pole1Tmp =   pole1StoreNew{lauf};
        % pole2Tmp = pole2StoreNew{lauf};
         
    redStackMid = I0File(:,:,lauf);
    redStackMid(:,:) = 0;
    
  cPoleMark1 =  cNewLongAxis1(lauf);
  rPoleMark1 =  rNewLongAxis1(lauf);
    
  cPoleMark2 = cNewLongAxis2(lauf);
  rPoleMark2 = rNewLongAxis2(lauf);
    
    
   cMark1 = cNewMid1(lauf);
   rMark1 = rNewMid1(lauf);
    
    
    cMark2 = cNewMid2(lauf);
    rMark2= rNewMid2(lauf);
    
    
    
    
    
    
    cFurrowContour1 =  flank1Tmp(:,1);
    rFurrowContour1 =  flank1Tmp(:,2);
    
    cFurrowContour2 =  flank2Tmp(:,1);
    rFurrowContour2 =  flank2Tmp(:,2);
     %%%%%
    redStackPole1 = I0File(:,:,lauf);
    redStackPole2 = I0File(:,:,lauf);
    
    
%     cPoleContour1 =    pole1Tmp(:,1);
%     rPoleContour1 =   pole1Tmp(:,2);
%     
%     
%     cPoleContour2 =    pole2Tmp(:,1);
%     rPoleContour2 =   pole2Tmp(:,2);
    %rMark1 = rMark1-2;
    %cMark1 = cMark1-2;
         redStackMid(rMark1,cMark1) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ cFurrowContour1';rFurrowContour1'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid1{lauf} = MIJ.getColumn('y');
         xRedMid1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redMax1(lauf) = find(yRedMid1{lauf}==max(yRedMid1{lauf}))
         
         redStackMid(:,:) = 0;
         redStackMid(rMark2,cMark2) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ cFurrowContour2';  rFurrowContour2'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid2{lauf} = MIJ.getColumn('y');
         xRedMid2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redMax2(lauf) = find(yRedMid2{lauf}==max(yRedMid2{lauf}))
         
%          
%             
%          %%%%%%% Max pole region
%         % cPoleContour1 = cPoleContour1-2;
%          %rPoleContour1 = rPoleContour1-2;
%          redStackPole1(:,:) = 0;
%          redStackPole1(rPoleMark1,cPoleMark1) = 65000;
%          MIJ.createImage(redStackPole1(:,:));
%          MIJ.run('setLine12');
%            MIJ.setRoi( [ cPoleContour1'; rPoleContour1'], ij.gui.Roi.POLYLINE);
%          MIJ.run('getLinescanRed');
%          yPoleMid1{lauf} = MIJ.getColumn('y');
%          xPoleMid1{lauf} = MIJ.getColumn('x');
%          MIJ.run('closeAllWindows');
%          
%          redPoleMax1(lauf) = find(yPoleMid1{lauf}==max(yPoleMid1{lauf}))
%          
%          redStackPole2(:,:) = 0;
%          redStackPole2(rPoleMark2,cPoleMark2) = 65000;
%          MIJ.createImage(redStackPole2(:,:));
%          MIJ.run('setLine12');
%          MIJ.setRoi( [ cPoleContour2'; rPoleContour2'], ij.gui.Roi.POLYLINE);
%          MIJ.run('getLinescanRed');
%          yPoleMid2{lauf} = MIJ.getColumn('y');
%          xPoleMid2{lauf} = MIJ.getColumn('x');
%          MIJ.run('closeAllWindows');
%          
%          redPoleMax2(lauf) = find(yPoleMid2{lauf}==max(yPoleMid2{lauf}))
%          
         
clear redStackMid   
end





%%%% set sliding window size
windowSize = 1; %%% in microns

sWSet = round(windowSize / voxelX_mum);

%%%%%%%%%  test whether the sliding window is odd or even

if mod(sWSet,2) == 1
    
    sWSet = sWSet +1;
else
    
end


%sWSet = 40;
sW = sWSet;
sW1 = sW; %%% in pixels
sW2 = sW;



clear  greenValuesSliding_1
clear  greenValuesSliding_2
clear  redValuesSliding_1
clear  redValuesSliding_2

 clear yI45_Sliding  yI90_Sliding yI135_Sliding yIB_Sliding

for lauf = 1:AnalysisEnd    

    yI0_1_tmp = yI0_1{lauf};
    yI45_1_tmp = yI45_1{lauf};
    yI90_1_tmp = yI90_1{lauf};
    yI135_1_tmp = yI135_1{lauf};
    yIB_1_tmp = yIB_1{lauf};
    
        %%%%% This section ensures that the sliding windows fits the
        %%%%% selection. If not the selection is cropped to fit the window
    
    %    redMax1(dist1 == 1) = 2;
        
        if redMax1(lauf)+ (sW1/2) > length(yI0_1_tmp)
            
            sW1 = round(length(yI0_1_tmp)-(redMax1(lauf)))-1
            if  mod(sW1,2) ==1
                sW1 = sW1 -1;
            end
            
        elseif redMax1(lauf)-(sW1/2) < 1
            sW1 = round(redMax1(lauf)/2)-1
            if  mod(sW1,2) ==1
                sW1 = sW1 -1;
            end
        end
        
   yI0_1_Sliding_1{lauf} = yI0_1_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yI45_1_Sliding_1{lauf} = yI45_1_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yI90_1_Sliding_1{lauf} = yI90_1_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yI135_1_Sliding_1{lauf} = yI135_1_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yIB_1_Sliding_1{lauf} = yIB_1_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
    

    yI0_2_tmp = yI0_2{lauf};
    yI45_2_tmp = yI45_2{lauf};
    yI90_2_tmp = yI90_2{lauf};
    yI135_2_tmp = yI135_2{lauf};
    yIB_2_tmp = yIB_2{lauf};
    
     if redMax2(lauf)+ (sW2/2) > length(yI0_2_tmp)
            
            sW2 = round(length( yI0_2_tmp)-(redMax2(lauf)))-1
            if  mod(sW2,2) ==1
                sW2 = sW2 -1;
            end
        elseif redMax2(lauf)-(sW2/2) < 1
             sW2 = round(redMax2(lauf)/2)-1
             if  mod(sW2,2) ==1
                sW2 = sW2 -1;
            end
        end
    sWStore1(lauf) = sW1
      sWStore2(lauf) = sW2
   
   yI0_2_Sliding_1{lauf} = yI0_2_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yI45_2_Sliding_1{lauf} = yI45_2_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yI90_2_Sliding_1{lauf} = yI90_2_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yI135_2_Sliding_1{lauf} = yI135_2_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
   yIB_2_Sliding_1{lauf} = yIB_2_tmp([redMax1(lauf)-(sW1/2):redMax1(lauf)+(sW1/2)]);
    
    
   
    
    
    %%%%%%%%%%%%%%%%%%
    
    
   
   
    %%%%%% Mean intensity per sliding windows
    
    
   yI0_Sliding(lauf) = mean((yI0_2_Sliding_1{lauf} + yI0_2_Sliding_1{lauf})/2)
   yI45_Sliding(lauf) = mean((yI45_2_Sliding_1{lauf} + yI45_2_Sliding_1{lauf})/2)
   yI90_Sliding(lauf) = mean((yI90_2_Sliding_1{lauf} + yI90_2_Sliding_1{lauf})/2)
   yI135_Sliding(lauf) = mean((yI135_2_Sliding_1{lauf} + yI135_2_Sliding_1{lauf})/2)
   yIB_Sliding(lauf) =mean((yIB_2_Sliding_1{lauf} + yIB_2_Sliding_1{lauf})/2)
    
    
    
    sW = sWSet;
sW1 = sW; %%% in pixels
sW2 = sW;

end



for lauf = 1:p
    ItoSMatrix=[0.5 0.5 0.5 0.5; 1 0 -1 0; 0 1 0 -1];
    S=ItoSMatrix*[yI0_Sliding(lauf)  yI45_Sliding(lauf)  yI90_Sliding(lauf) yI135_Sliding(lauf)]';

    aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
         aniso=(1+aniso)./(1-aniso); %Then Ratio
         aniso(aniso<1)=NaN;
         
         ratioStorePole(lauf) = aniso
         
end



for lauf = 1:p
  S=ItoSMatrix*[yI0_Sliding(lauf)  yI45_Sliding(lauf)  yI90_Sliding(lauf) yI135_Sliding(lauf)]';
     aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);

   
         
         anisoStore(lauf) = aniso
         
end



h= figure
plot(ratioStorePole,'r')
hold on
plot(anisoStore,'b')

%axis([timeVec(1) 1000 0 25])  
%yL = get(gca,'YLim');
%line([ingressionTime ingressionTime],yL,'Color','r');
 
legend('Ratio cleavage furrow','Ratio pole')

 
xlabel ('Time [s]','FontSize', 16);
ylabel('Ratio [A.U.]','FontSize', 16);
title('Polarization ratio polscope','FontSize', 16);

print(h,'-dpdf', ['_polscopeRatio.pdf']);%tifCurvetifFilename);

                
                
                for lauf = 1:p
                    
                    
                    
                    [BW_Box BWI0 BWI45 BWI90 BWI135 BWI135B x y...
                        I0Avg I45Avg I90Avg I135Avg I0BleachAvg] = doDrawAnalysisBox5Orientations(I0File(:,:,lauf),I0File(:,:,lauf), I45File(:,:,lauf),...
                        I90File(:,:,lauf), I135File(:,:,lauf), I0File(:,:,lauf) ,voxelX_mum)
                    
                    
                    
                    xFurrow1{lauf} = x;
                    yFurrow1{lauf} = y;
                    
                    
                    
                    
                    %%%% Parse data
                    
                    I0_Furrow(lauf) =(I0Avg(1) + I0Avg(2))/2;
                    I0_Pole(lauf) = (I0Avg(3) + I0Avg(4))/2;
                    
                    I45_Furrow(lauf) = (I45Avg(1) + I45Avg(2))/2;
                    I45_Pole(lauf) = (I45Avg(3) + I45Avg(4))/2;
                    
                    I90_Furrow(lauf) = (I90Avg(1) + I90Avg(2))/2;
                    I90_Pole(lauf) = (I90Avg(3) + I45Avg(4))/2;
                    
                    I135_Furrow(lauf) = (I135Avg(1) + I135Avg(2))/2;
                    I135_Pole(lauf) = (I135Avg(3) + I135Avg(4))/2;
                    
                    I0Bleach_Furrow(lauf) = (I0BleachAvg(1) + I0BleachAvg(2))/2;
                    I0Bleach_Pole(lauf) = (I0BleachAvg(3) + I0BleachAvg(4))/2;
                    
                    
                    I0_Store{lauf} = I0Avg;
                    I45_Store{lauf} = I45Avg;
                    I90_Store{lauf} = I90Avg;
                    I135_Store{lauf} = I135Avg;
                    I0Bleach_Store{lauf} = I0BleachAvg;
                    BW_BoxStore{lauf} =    BW_Box;
                    
                    BWI0_Store{lauf} = BWI0;
                    BWI45_Store{lauf} = BWI45;
                    BWI90_Store{lauf} = BWI90;
                    BWI135_Store{lauf} = BWI135;
                    BWI135B_Store{lauf} = BWI135B;
                    x_Store{lauf} = x;
                    y_Store{lauf} = y;
                end

          for lauf = 1:p
              x_StoreTmp = cell2mat(x_Store{lauf});
              y_StoreTmp = cell2mat(y_Store{lauf});
              
               cleavageFurrowDiameter(lauf) = (pdist2([x_StoreTmp(1),y_StoreTmp(1)],[x_StoreTmp(2),y_StoreTmp(2)])) * voxelX_mum
          end
          
          plot(cleavageFurrowDiameter)

      
          

plot(I0_Furrow)
            



arg.Zslices=1; 
arg.Frames=1;
arg.Channel=1; % Select the channel to process.
%arg.Position=1;
arg.registerPol=false; % Do pol-channels need registration?
arg.bin=1;
arg.normFactors=[1 1 1 1];
arg.suffix=''; %Suffix after the slice name.
arg.prefix=''; %Prefix before the time-stamp for data acquired by polling when files are written by Zen.
arg.BlackLevel=0;
arg.bitDepth=NaN; %Bit-depth has no effect on how data is read. It only affects the format in which computed results are written.
arg.displayStatus=true;
arg.acqMethod='zen-controller';
%arg=parsepropval(arg,varargin{:});
clear anisoStore
clear avgStore
for lauf = 1:AnalysisEnd
    
  %  [orient, aniso, avg]=...
  %      ComputeFluorAnisotropy(I0_Furrow(lauf),I45_Furrow(lauf),I90_Furrow(lauf),I135_Furrow(lauf),...
  %      'anisotropy','BlackLevel',arg.BlackLevel,'normFactors',arg.normFactors);
    [orientCorr, anisoCorr, avgCorr]=ComputeFluorAnisotropy(I0_Furrow(lauf),I45_Furrow(lauf),I90_Furrow(lauf),I135_Furrow(lauf),...
        'anisotropy','I0bleach',I0_Furrow(lauf),'bleachROI',1,'normFactors',normFactors);

    %orientStore{lauf} = orient;
    anisoStore(lauf) = anisoCorr;
    avgStore(lauf) = avgCorr;
    
end

for lauf = 1:p
    ItoSMatrix=[0.5 0.5 0.5 0.5; 1 0 -1 0; 0 1 0 -1];
    S=ItoSMatrix*[I0_Furrow(lauf) I45_Furrow(lauf) I90_Furrow(lauf) I135_Furrow(lauf)]';

    aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
         aniso=(1+aniso)./(1-aniso); %Then Ratio
         aniso(aniso<1)=NaN;
         
         ratioStore(lauf) = aniso
         
end


for lauf = 1:p
    ItoSMatrix=[0.5 0.5 0.5 0.5; 1 0 -1 0; 0 1 0 -1];
    S=ItoSMatrix*[I0_Pole(lauf) I45_Pole(lauf) I90_Pole(lauf) I135_Pole(lauf)]';

    aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
         aniso=(1+aniso)./(1-aniso); %Then Ratio
         aniso(aniso<1)=NaN;
         
         ratioStorePole(lauf) = aniso
         
end



for lauf = 1:p
  S=ItoSMatrix*[I0_Furrow(lauf) I45_Furrow(lauf) I90_Furrow(lauf) I135_Furrow(lauf)]';
     aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);

   
         
         anisoStore(lauf) = aniso
         
end



h= figure
plot(ratioStore,'r')
hold on
plot(ratioStorePole,'b')

%axis([timeVec(1) 1000 0 25])  
%yL = get(gca,'YLim');
%line([ingressionTime ingressionTime],yL,'Color','r');
 
legend('Ratio cleavage furrow','Ratio pole')

 
xlabel ('Time [s]','FontSize', 16);
ylabel('Ratio [A.U.]','FontSize', 16);
title('Polarization ratio polscope','FontSize', 16);

print(h,'-dpdf', ['_polscopeRatio.pdf']);%tifCurvetifFilename);



for lauf = 1:AnalysisEnd
    
[orientCorr, anisoCorr, avgCorr]=ComputeFluorAnisotropy(I0,I45b,I90b,I135b,'anisotropy','I0bleach',I0b,'bleachROI',1,'normFactors',normFactors);
   % [orient, aniso, avg]=...
    %    ComputeFluorAnisotropy(I0_Pole(lauf),I45_Pole(lauf),I90_Pole(lauf),I135_Pole(lauf),...
    %    'ratio','BlackLevel',arg.BlackLevel,'normFactors',arg.normFactors);
    
    %orientStore{lauf} = orient;
    anisoStorePole(lauf) = anisoCorr;
    avgStorePole(lauf) = avgCorr;
    
end

ratioFurrowPole = ratioStore./ratioStorePole

figure(1)
figure(2)
plot(ratioFurrowPole)
hold on

plot(ratioStore)

ratioStore = anisoStore
%
close all
plot(anisoStorePole)
plot(avgStorePole)
plot(avgStore)


plot(avgAverage)

[AX,H1,H2] = plotyy(1:AnalysisEnd,anisoStore,1:AnalysisEnd,avgStore);

[AX,H1,H2] = plotyy(1:AnalysisEnd,anisoStore,1:AnalysisEnd,cleavageFurrowDiameter);

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
    
        
   
   cd(curdir)
   
       save([saveFileName,'.mat'])
   
            dummyVariable = (1:p)';
            dummyVariable(:) = 0;
           

            saveVariables = {};
            saveVariables{1} = timeVec';
            saveVariables{2} = cleavageFurrowDiameter';
            saveVariables{3} =  avgStore';
            
            saveVariables{4} =   avgStorePole';
             saveVariables{5} =  anisoStore';
            saveVariables{6} =   anisoStore';
            
           
            
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, cleavageFurrowDiameter, intensityCleavageFurrow, intensityPole',...
                ',ansioCleavageFurrow, anisoPole']
            
            outid = fopen([saveFileName,'Analysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([saveFileName,'Analysis.csv'],csvData,'roffset',1,'-append')
            
            
   
   
   
   
   
   
   
   
   [anisoPath,orientPath,avgPath]=writeComputedChannels(aniso,orient,avg,directory,idZ,idF,arg.Channel,arg.suffix,bitDepth);
        filenames=cat(1,filenames,{anisoPath,orientPath,avgPath}');

  
        
        
       