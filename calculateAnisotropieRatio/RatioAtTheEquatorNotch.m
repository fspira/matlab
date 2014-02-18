%%%%%% CaluclateAnisotropie at the cell equator. This script will first
%%%%%% segment the cell - then allow the user to validate segmentation.
%%%%%% Second user selection of centrall furrow region as well as a region
%%%%%% with simillar curvatore outside the central furrow. In addition
%%%%%% chromatin-chromatin distance hast to be selected by the user Third, curvature
%%%%%% is computed and furrow vs. non furrow curvature is computed - to
%%%%%% assure comparable curvature. Fourth, computed furrow ingression
%%%%%% diameter and chromatin-chromatin distance.


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
addpath('/Users/spira/Desktop/programme')
addpath('/Users/spira/Documents/Matlab_scripte/tiffIO')
addpath('/Users/spira/Documents/Matlab_scripte/')
addpath('/Users/spira/Desktop/Desktop/LifeactCherry_GlGPIEgfp/131204')
addpath('/Users/spira/Documents/MATLAB_scripte/ImageProcessing/Utilities')
curdir = pwd;

%load('voxelX_mum.mat');
%load('voxelX_mumMid');
 
tifFilename = 'cell12_FullStack.lsm';
tifFilenameMid = 'cell12_mid.lsm';

zSectionToAnalyze = 19% notch casette image
zSectionMidStack = 1; % conventional image

%%%%% Load midSection

imgMid = imread(char(tifFilenameMid));
imgMidtmp = tiffread30(char(tifFilenameMid));

imgMidtmpTmp = cat(3,imgMidtmp.data);
imgMid = imgMidtmpTmp;
imgMidtmpTmp = imgMidtmpTmp(:,:,zSectionMidStack);

voxelX = getfield(imgMidtmp,'lsm','VoxelSizeX');
voxelX_mumMid = voxelX*1000000;

%for lauf = 1:3
    
%    imgMid(:,:,lauf) = imgMidtmpTmp{lauf};
    
%end

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
imgOrigGreen = img{1};
imgOrigRed = img{2};

[m n p] = size(img);

p=1

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

[redImg greenImg] = bkgCorrectionRedGreen(redImg, greenImg);

%%%%% Setting negative values to small values, to avoid division through zero
greenNormZero = greenImg > 0;
greenNorm = double((greenImg .* greenNormZero));

redNormZero = redImg > 0;
redNorm = double(redImg .* redNormZero);

%%%%% change directory
cd ([curdir '/' folderName]);
workdir = pwd;


S1 = greenNorm;
S2 = redNorm;

D=((S1-1.079.*S2)./(S1+(2*1.079).*S2)).*255;
D(~isfinite(D)) = 0;
[D, dNegInv, dSum] = doVisualizeSPComponents(D,tifFilename);

%%%%%% normalize to mean itensity of each channel
[m n p] = size(img);
for lauf = 1:p
    S1Norm(:,:,lauf) = S1(:,:,lauf)./mean(mean(S1(:,:,lauf)));
    S2Norm(:,:,lauf) = S2(:,:,lauf)./mean(mean(S2(:,:,lauf)));
    ratio(:,:,lauf) = S2Norm(:,:,lauf)./S1Norm(:,:,lauf);
  
    
end


segmentOut = doSegmentImage(S1);

   
  cc = regionprops(segmentOut,'all');
  B = bwboundaries(segmentOut,8); 
 B  = B{1};


%%%%% manually confirm segmentation


Miji;

cd(curdir)


                MIJ.createImage(imgMerge(:,:,:,lauf));
                MIJ.run('8-bit');
              
               MIJ.run('Stack to RGB');
               % MIJ.selectWindow('Import from Matlab');
              MIJ.setRoi( [B(:,2)';B(:,1)'], ij.gui.Roi.POLYLINE);
      %MIJ.setRoi( [ selectedRoi(:,2)'; selectedRoi(:,1)'], ij.gui.Roi.POLYLINE);
    k = waitforbuttonpress 
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    cSegment{lauf} = tmp(:,1);
    rSegment{lauf} = tmp(:,2);
    
    segmentedOutline = [cSegment{1}, rSegment{1}];
    
    BW1 = S2;
    BW1(:,:) = 0;
    
   for subrun = 1:length(rSegment{lauf}) 

    BW1(segmentedOutline(subrun,2), segmentedOutline(subrun,1)) = 255;
    
   end
   
   
   
   imshow(BW1)
   
   [splineFitOutline, imgOut ]= doASplineFit(segmentedOutline,S1);

   ccc = regionprops(imgOut,'Orientation','MajorAxisLength','MinorAxisLength','Centroid');

   k=1
   orientHor(1) = 90;
   orientHor(2) = -90;
   orientVer(1) = 0;
   orientVer(2) = 180;
   
   cellOrientation = sqrt(ccc.Orientation^2)
   
       if cellOrientation > 45

           for subrun = 1:2
                
               [xp(subrun) yp(subrun)] = doSelectEllipsePoint(ccc, orientHor(subrun) ,k)    
                
           end
           
       else
           
           for subrun = 1:2
                
                [xp(subrun) yp(subrun)] = doSelectEllipsePoint(ccc, orientVer(subrun) ,k)    
                
           end
           
       end
       
       
       %%%%%%% doDrawCenterOfFurrow
       
         

       [Furrow1 Furrow2 Curve1 Curve2]= doPickFurrowSecondCurve(imgMerge)
               
          dist1 = pdist2([Furrow1(1),Furrow1(1)],[Furrow2(1),Furrow2(2)]);
                   
         furrowDiameter = dist1*voxelX_mum;
       
          [furrowCenter1,furrowCenter2, curveCenter1, curveCenter2]= doDetectCenterOfRoi(redNorm ,splineFitOutline,Furrow1,Furrow2,Curve1,Curve2)
          [S1Linescan, S2Linescan] = doGetLinescan(S1Norm, S2Norm,splineFitOutline);
          
          S1Linescan = S1Linescan{1};
          S2Linescan = S2Linescan{2};
          
        %%%% set sliding window size in micron
windowSize = 4; %%% in microns

sWSet = round(windowSize / voxelX_mum);  

  furrow1Sliding_S1 = getIntensityUnderSlidingWindow(S1Linescan,sWSet,furrowCenter1)
  furrow1Sliding_S2 = getIntensityUnderSlidingWindow(S2Linescan,sWSet,furrowCenter1)
  furrow2Sliding_S1 = getIntensityUnderSlidingWindow(S1Linescan,sWSet,furrowCenter2)
  furrow2Sliding_S2 = getIntensityUnderSlidingWindow(S2Linescan,sWSet,furrowCenter2)
       
  
  curve1Sliding_S1 = getIntensityUnderSlidingWindow(S1Linescan,sWSet,curveCenter1)
  curve1Sliding_S2 = getIntensityUnderSlidingWindow(S2Linescan,sWSet,curveCenter1)
  curve2Sliding_S1 = getIntensityUnderSlidingWindow(S1Linescan,sWSet,curveCenter2)
  curve2Sliding_S2 = getIntensityUnderSlidingWindow(S2Linescan,sWSet,curveCenter2)
  
  
  
  ratioScanS1S2 = S1Linescan ./ S2Linescan;
  
  ratioSliding_Furrow1 = getIntensityUnderSlidingWindow(ratioScanS1S2,sWSet,furrowCenter1)
  ratioSliding_Furrow2 = getIntensityUnderSlidingWindow(ratioScanS1S2,sWSet,furrowCenter1)
  
  plot(1:length(ratioSliding_Furrow1), ratioSliding_Furrow1,'-')
  yL = get(gca,'YLim');
  hold on
  line([furrowCenter1 furrowCenter1],yL,'Color','r');
  line([furrowCenter2 furrowCenter2],yL,'Color','b');
  line([curveCenter1 curveCenter1],yL,'Color','g');
  line([curveCenter2 curveCenter2],yL,'Color','g');

    
    
  plot(1:length( ratioScanS1S2),  ratioScanS1S2,'-')
  yL = get(gca,'YLim');
  hold on
  line([furrowCenter1 furrowCenter1],yL,'Color','r');
  line([furrowCenter2 furrowCenter2],yL,'Color','b');
  line([curveCenter1 curveCenter1],yL,'Color','g');
  line([curveCenter2 curveCenter2],yL,'Color','g');
  
  
  % connectedPoints = doConnectTwoPoints(xp,yp,imgOut)
    %   yCat = doInterpolateLineThroughPoints(xp, yp, imgOut)
       
    %   tmp1 = imgOut; tmp1(:,:)=0; %tmp1 = im2bw(tmp1);
       
    %   for subrun = length(yCat)
          
    %      tmp1((round(yCat(subrun,1))),round(yCat(subrun,2))) = 125;
            
           
    %   end
           
          %  tmp1 = bwmorph(tmp1,'bridge',8);
          %  tmp1 = im2bw(tmp1);
          %  tmp1 = tmp1.*125;
         
   
    

    %tmp2 = tmp + tmp1;
    %imshow(tmp2)


%[rIntersect2,cIntersect2] = find(tmp2 > 135)
           

