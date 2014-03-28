%%%%%% Pogram to draw a kymograph over a movie. This program will segment
%%%%%% the image and ask the user to update segmentation result. The
%%%%%% linescaen will be used for the kymograph
%%%%%%
%%%%%%


clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

blue = zeros(256,3);
blue(:,3) = distvec;

javaaddpath '/Applications/Fiji.app/scripts/Miji.m'
javaaddpath '/Applications/Fiji.app/jars/ij-1.48s.jar'

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/mij.jar'


addpath('/Applications/Fiji.app/scripts')

addpath('/Applications/Fiji.app/plugins')
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
Miji;
%MIJ.start('/Applications/Fiji')

anaOnset = 3;
ingressionFrame = 9;
lastFrameToConsider = 16;

curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell8_100nm_vertical.lsm';
%tifFilenameMid = 'cell6_conv.tif';

saveFileName = 'cell8_100nm_vertical_center_CurvatureKymo';

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

%%%%%%
%%%%%% Segment the image
%%%%%%

    boundary = doSegmentImage(redImg);

    %%%%%
    %%%%% validate the segmentaiotn
    %%%%%
   
    
for lauf = 1: lastFrameToConsider
    
     L2 = boundary(:,:,lauf);
  cc = regionprops(L2,'all');
  B = bwboundaries(L2,8); 
  boundaryCoords = B{1};
 
    MIJ.createImage(redImg(:,:,lauf));
    MIJ.run('8-bit');
  %  MIJ.run('Stack to RGB');
   % MIJ.selectWindow('Import from Matlab');
    MIJ.setRoi( [boundaryCoords(:,2)';boundaryCoords(:,1)'], ij.gui.Roi.POLYLINE);
      %MIJ.setRoi( [ selectedRoi(:,2)'; selectedRoi(:,1)'], ij.gui.Roi.POLYLINE);
    k = waitforbuttonpress 
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    cNew1Ref{lauf} = tmp(:,1);
    rNew1Ref{lauf} = tmp(:,2);
    
   MIJ.run('closeAllWindows');
end

%%%%%% 
%%%%%% generate ROI for each flank of the furrwo
%%%%%%
   
channelMergeImage = (redImg(:,:,:)+ greenImg(:,:,:))./2;
    
for lauf = 1:lastFrameToConsider
     [xCoord,yCoord,regionIntersect,furrowRoi1,LinescanFurrow1]= splitFurrows(redImg(:,:,lauf),xFurrow1{lauf},yFurrow1{lauf},cNew1Ref{lauf},rNew1Ref{lauf},150,channelMergeImage(:,:,lauf)) ;   

     
     LinescanFurrow1_Store{lauf} = LinescanFurrow1;
     furrow1MidCoords{lauf} = [xCoord yCoord];
     furrowRoi1_Sammel{lauf} = furrowRoi1;
     midFurrow1(lauf) = regionIntersect;

      [xCoord,yCoord,regionIntersect,furrowRoi2,LinescanFurrow2]= splitFurrows(redImg(:,:,lauf),xFurrow2{lauf},yFurrow2{lauf},cNew1Ref{lauf},rNew1Ref{lauf},150,channelMergeImage(:,:,lauf)) ;
     
      LinescanFurrow2_Store{lauf} = LinescanFurrow2;
      furrowRoi2_Sammel{lauf} = furrowRoi2;
     
        furrow2MidCoords{lauf} = [xCoord yCoord];
      midFurrow2(lauf) =  regionIntersect;
end


%%%%%%%
%%%%%%% kymograph section
%%%%%%%

%%%%% Cut linescan 140 px around the center pixeld



figure(1)
hold on
 for lauf = 1:lastFrameToConsider
    centerPixel_1 =  midFurrow1(lauf)
    LinescanFurrow1Tmp_1 = LinescanFurrow1_Store{lauf}
    LinescanFurrow1LinearTmp_1 = LinescanFurrow1Tmp_1(centerPixel_1-70:centerPixel_1+70);
    
    centerPixel_2 =  midFurrow2(lauf)
    LinescanFurrow1Tmp_2 = LinescanFurrow2_Store{lauf}
    LinescanFurrow1LinearTmp_2 = LinescanFurrow1Tmp_2(centerPixel_2-70:centerPixel_2+70);
    LinescanMerge = (LinescanFurrow1LinearTmp_2 + LinescanFurrow1LinearTmp_1)./2;
    
    LinescanLinearStore(lauf,:) = LinescanMerge';
   
    Linescan1(lauf,:) = LinescanFurrow1LinearTmp_1';
    Linescan2(lauf,:) = LinescanFurrow1LinearTmp_2';
   
 end
 
  for lauf = 1:lastFrameToConsider
    plot(1:length(LinescanLinearStore(lauf,:)),LinescanLinearStore(lauf,:),'-')
    hold on
   % plot(centerPixel_2, LinescanFurrow1LinearTmp_2(centerPixel_2),'xr')
    
    hold off
    pause(0.5)
    
  end
  
  
  %%%%%%% Determine the analysis window
  %%%%%%%
  %%%%%%% use the mean of the first and last values of the ROI - this
  %%%%%%% should represent values not part of the accumulation area. Use
  %%%%%%% intensities higher than BKG for analysis.
  
  for lauf = 1:lastFrameToConsider
      
   
    LinescanTmp =  LinescanLinearStore(lauf,:);
    LinescanBkg = (mean(LinescanTmp(1:10))+mean(LinescanTmp(end-10:end)))/2;
       centerPixel =round(length(LinescanTmp)/2);
    %%%%%% Add 10% to Bkg to ensure that the signal chosen for analysis is
    %%%%%% real
    LinescanBkgSub = LinescanTmp - (LinescanBkg+((LinescanBkg/100)*10));
    
    plot(1:length(LinescanBkgSub), LinescanBkgSub)
    %%%%%% Eliminate elements smaller than zero
    LinescanBkgSubTmp = LinescanBkgSub >= 0;
   
    runMarker = 1;
    clear leftEdge;
    clear rightEdge;
    
    for subrun = 1:centerPixel-1
        
      if  LinescanBkgSubTmp(centerPixel-subrun) == 1
        
      else
            leftEdge(runMarker) = centerPixel-subrun
            runMarker = runMarker +1
      end
            
    end
      runMarker = 1;
     for subrun = 1:  (length(LinescanBkgSubTmp) - centerPixel)-1
        
       if LinescanBkgSubTmp(centerPixel+subrun) == 1
        
       else
            rightEdge(runMarker) = centerPixel+subrun
             runMarker = runMarker +1;
       end
            
     end
    
     LinescanAnalysis = LinescanBkgSub( leftEdge(1):rightEdge(1));
    
     plot(1:length(LinescanAnalysis), LinescanAnalysis)
    
    
  if length(LinescanAnalysis) <= 20
      
      FWHM_Store{lauf} = 0;
      y_fitted_Store{lauf} = 0;
      x_ax_Store{lauf} = 0;
      maxDistance_Store{lauf} = 0;
      
  else

           distance = 0+voxelX_mum:voxelX_mum: length(LinescanAnalysis) * voxelX_mum;
           maxDistance = length(LinescanAnalysis) * voxelX_mum;

          [fittedmodel, goodness, output] = fit(distance',LinescanAnalysis','gauss1');
          sigma1 = fittedmodel.c1/sqrt(2); %%%% Gives sigma
          FWHM = 2*sqrt(2*log(2))*sigma1 %%% Calculates the FWHM

        x_ax = -400:0.01:400; % axis
        y_fitted = fittedmodel(x_ax); 


          fwhmTest1 = fwhm(x_ax,y_fitted)

        h = figure

        plot(distance,LinescanAnalysis,'r');
        hold on
        plot(x_ax,y_fitted)
        axis([-4 14 0 max(LinescanAnalysis)+40])  


        title(['FWHM', char(lauf)])
        xlabel ('Distance [µm]','FontSize', 16);
         ylabel('Intensities [A.U.]','FontSize', 16);



        print(h,'-dpdf', [curdir '/' tifFilename,'_FWHM.pdf']);%tifCurvetifFilename);
        close all

         
      FWHM_Store{lauf} = FWHM;
      y_fitted_Store{lauf} = y_fitted;
      x_ax_Store{lauf} = x_ax;
      maxDistance_Store{lauf} = maxDistance;
  
 % provides all the fitted values for  
                        %%%%the specified x-axis

    end
  end

  for lauf = 1:lastFrameToConsider
      plot(lauf,maxDistance_Store{:,lauf},'bx')
      hold on
       plot(lauf,FWHM_Store{:,lauf},'rx')
  end
 
 %%%%% This piece of code determines the axis of the cell and suggests
    %%%%% position of pole-pole and furrow-furrow centers
   % lauf = 15
   
    img= greenImg(:,:,1);
   img(:,:) = 0;
 for   lauf = 1:lastFrameToConsider;
   newFurrow = furrowRoi2_Sammel{lauf}
  MIJ.createImage(channelMergeImage(:,:,1));
         %MIJ.run('setLine8');
     % MIJ.setRoi( [ boundaryRoi(:,1)';  boundaryRoi(:,2)'], ij.gui.Roi.POLYLINE);
      MIJ.setRoi( [ newFurrow(:,1)'; newFurrow(:,2)'], ij.gui.Roi.POLYLINE);
       MIJ.run('getLinescanRed');
      yRedAxisTest = MIJ.getColumn('y')
            MIJ.closeAllWindows

 end
 
 MIJ.createImage( LinescanLinearStore)
 

%%%%
 [Curvature_Furrow1] = doCalculateCurvatureAnisotropie(furrowRoi1_Sammel, 80, midFurrow1,greenImg(:,:,1),furrow1MidCoords)

 [Curvature_Furrow2] = doCalculateCurvatureAnisotropie(furrowRoi2_Sammel, 80, midFurrow2,greenImg(:,:,1),furrow2MidCoords)
 
 CurvatureMean = (Curvature_Furrow1(:,3)+Curvature_Furrow2(:,3))./2;


plot(1:lastFrameToConsider,Curvature_Furrow2(:,3),'rx')
hold on
plot(1:lastFrameToConsider,Curvature_Furrow1(:,3),'bx')
plot(1:lastFrameToConsider,CurvatureMean,'cx')

 axis([0 lastFrameToConsider 0 200])   

hold off

plot(1:lastFrameToConsider,CurvatureMean,'rx')
 axis([0 lastFrameToConsider 0 500])   
 
 for lauf =1:lastFrameToConsider
     imshow(greenImg(:,:,lauf),[])
     hold on
     circle(Curvature_Furrow1(lauf,1),Curvature_Furrow1(lauf,2),Curvature_Furrow1(lauf,3))
     circle(Curvature_Furrow2(lauf,1),Curvature_Furrow2(lauf,2),Curvature_Furrow2(lauf,3))
     pause(0.5)
  hold off
 end
 
 
 
 
 
 
