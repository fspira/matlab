%%%%% This script uses the workspace of intensity measuremnts as an input.
%%%%% It requires the linescan of two channels (can be adapted for single
%%%%% channel usage) and the center position. The script will normalize the
%%%%% curves by subtracting the maean of 5% of left and right from the
%%%%% linescan. Then it will fit a single gaussian to the center position
%%%%% (10pixel variation allowed as starting value) and calculates the
%%%%% FWHM.



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
addpath('/Users/spira/Desktop/programme/tools/bf')


FWHMStore ={};
FWHMOrigStoreGreen =[];
 FWHMOrigStoreRED = [];
 FWHMTimeVec = [];
 FWHMOrigStore = [];

%%%% Average flanks
[greenOut redOut midOut] = doAverageFlanks(redMax1,redMax2,yGreen1,yGreen2,yRed1,yRed2)


curdir = pwd;
%%%% Select frame to plot

%%%% Select the center position
maxSelector = midOut;

%%%% Select the linescanchannel
linescanSelectorIn = greenOut;


%%%% Use the poles as a reference for normalization
%%%% Use the total intensity for normalization

%for lauf = 1:length(linescanSelectorIn)

 %   linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yGreenNonPol1{lauf})+mean(yGreenNonPol2{lauf}) + mean(yGreen1{lauf}) + mean(yGreen2{lauf}))/4;

%end


%%% Smooth the curve a sliding average of 5 is used
%%% Since the window creates artifacts at the left edge, a dummy value is
%%% first attached to the left edge and (mean of 5 Values and subsequently
%%% removed

%windowSize = 5;

%b = (1/windowSize)*ones(1,windowSize)
%a = 1;

%for lauf = 1:length(linescanSelectorIn)
%    linescanTmp = linescanSelectorIn{lauf};
%    linescanTmpMean = mean(linescanTmp(1:5))
%    linescanFill(1:5) = linescanTmpMean ;
%    linescanTmp = cat(1,linescanFill',linescanTmp);
    
 %   linescanSelectorTmp = filter(b,a,linescanTmp);
 %   linescanSelectorIn{lauf} = linescanSelectorTmp(6:length(linescanSelectorTmp))
    % plot(length(linescanSelectorIn{lauf}),linescanSelectorIn{lauf},'b')
%end

%%%%% use a flatline background correction - the left and right 20 pixels
%%%%% are used for correction

%for lauf = 1:length(linescanSelectorIn)
 %   linescanTmp = linescanSelectorIn{lauf};
 %   [y,yfit] = bf( linescanTmp,[1:5,length(linescanTmp)-5:length(linescanTmp)],'linear');
    
  %  y_Indx = y > 0;
  %  yNew = y_Indx .* y;

   % linescanSelectorIn{lauf} = yNew;
    
%end

%for lauf = 1:length(linescanSelectorIn)
%    plot(1:length(linescanSelectorIn{lauf}),linescanSelectorIn{lauf},'b')
   % plot(1:length(linescanSelector{lauf}),linescanSelector{lauf},'b')
      %axis([-20 +20 0 1]) 
 %     pause(0.2)
%end
%hold on
%figure(2)
%plot(testy,y,'r')

%%%% Detect minimum Values
%plot(1:length(linescanSelector{lauf}), linescanSelector{lauf})

%for lauf = 1:length(linescanSelectorIn)

  % linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yGreenNonPol1{lauf})+mean(yGreenNonPol2{lauf}) )/ 2;
   
   %linescanSelector{lauf} = linescanSelector{lauf} - min(linescanSelector{lauf}); 
   
   
%end


%for lauf = 1:length(linescanSelectorIn)

 %  linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yRedNonPol1{lauf})+mean(yRedNonPol2{lauf}) )/ 2;
  % linescanSelector{lauf} = linescanSelector{lauf} - min(linescanSelector{lauf}); 
   
   
%end

%%%%%% Normalize by subtracting a fixed value which is the average of 5 percent of the
%%%%%% lateral data points



for lauf = 1:length(linescanSelectorIn)

   linescanSelectorTmp =  linescanSelectorIn{lauf} ;
   averagePixel = round((length(linescanSelectorTmp)/100)*5);
    linescanAvgLeft = linescanSelectorTmp(1:averagePixel);
    linescanAvgRight =linescanSelectorTmp(end-averagePixel:end);
    
    linescanSubMean = mean((linescanAvgLeft)+mean(linescanAvgRight))/2;
    
  linescanSelector{lauf} = linescanSelectorIn{lauf} - linescanSubMean; 
   
end


%%%% Use the average of the total intensity for normalization
%%%% Has to be implemented
for lauf = 1:length(linescanSelectorIn)
 
lauf
    
    %%%% find frame relative to anaphase onset
    time1 = timeVec(lauf);
    
    
    
    %%%% load the data from the linescan and find the center position
    redMaxTmp = maxSelector(lauf)-maxSelector(lauf);
    yGreenTmp = linescanSelector{lauf}%- min(linescanSelector{lauf})
    %%%% Calculated the full width half maximum
    [FWHMOrig, x_axOrig, y_fittedOrig, centeredCurve] = doFWHM_Linescan( linescanSelector{lauf},maxSelector,voxelX_mum,lauf) ;

   
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

for lauf = 1: length(midOut)
    if midOut(lauf) < 150


        linescanSelectorTmp_I = linescanSelector{lauf}

       % length(linescanSelectorTmp_I);
       % midOut(lauf);

        addZero = zeros(1,150 - midOut(lauf))';

        linescanSelectorCrop = cat(1,addZero,linescanSelectorTmp_I);
        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorCrop,addZero);

        linescanSelectorCropStoreGreen(:,lauf) = linescanSelectorCrop(1:300);
    elseif midOut(lauf) > 150

         linescanSelectorTmp_I = linescanSelector{lauf}


        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreen(:,lauf) =  linescanSelectorCrop(midOut(lauf)-149:300+(midOut(lauf)-150));


    elseif midOut(lauf) == 150

        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1, linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreGreen(:,lauf) = linescanSelectorCrop(1:300);

    end
end

 FWHMOrigStoreGreen =  FWHMOrigStore';
 
 FWHMTimeVec= timeVec'
    %%%%% Fit cropped 10% above background
    %plot(x_ax,y_fitted,'r')

    %%%%% Fit non cropped curves
    %plot(x_axOrig,y_fittedOrig,'r')
     mkdir([curdir '\' 'redFit']);
       mkdir([curdir '\' 'greenFit']);

    for lauf = 20:length(linescanSelectorIn)
         redMaxTmp = maxSelector(lauf)-maxSelector(lauf);
        yGreenTmp = linescanSelector{lauf}
        h = figure(1) 
        plot(((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum),yGreenTmp,'c')
        hold on
       plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
          axis([-20 +20 0 500]) 
  title('Gaussian fit to MRLIIC width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16);
    ylabel('Intensities [A.U.]' ,'FontSize', 16);
    legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(FWHMTimeVec(lauf)) 's' ])
        pause(0.2)
       
       
      % print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5percAVG.pdf']);%tifCurvetifFilename);
   
       % print(h,'-dpdf', [curdir '\' 'redFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_noAverage_NopxBFminAvgSub.pdf']);%tifCurvetifFilename);
   close all
    end
    
   %%%%%%%%%%%%%%%%% RED CHANNEL


%%%%%% Normalize by subtracting a fixed value which is the average of 5 percent of the
%%%%%% lateral data points

linescanSelectorIn = redOut;

for lauf = 1:length(linescanSelectorIn)

   linescanSelectorTmp =  linescanSelectorIn{lauf} ;
   averagePixel = round((length(linescanSelectorTmp)/100)*5);
    linescanAvgLeft = linescanSelectorTmp(1:averagePixel);
    linescanAvgRight =linescanSelectorTmp(end-averagePixel:end);
    
    linescanSubMean = mean((linescanAvgLeft)+mean(linescanAvgRight))/2;
    
  linescanSelector{lauf} = linescanSelectorIn{lauf} - linescanSubMean; 
   
end


%%%% Use the average of the total intensity for normalization
%%%% Has to be implemented
for lauf = 1:length(linescanSelectorIn)
 
lauf
    
    %%%% find frame relative to anaphase onset
    time1 = timeVec(lauf);
    
    
    
    %%%% load the data from the linescan and find the center position
    redMaxTmp = maxSelector(lauf)-maxSelector(lauf);
    yGreenTmp = linescanSelector{lauf}%- min(linescanSelector{lauf})
    %%%% Calculated the full width half maximum
    [FWHMOrig, x_axOrig, y_fittedOrig, centeredCurve] = doFWHM_Linescan( linescanSelector{lauf},maxSelector,voxelX_mum,lauf) ;

   
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
    centeredCurveStoreRed{lauf} = [centeredCurve' linescanSelector{lauf}]
  %  plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
end


%%%%%% This part shifts the linescan and centeres the linescan at zero
%%%%%% position

for lauf = 1: length(midOut)
    if midOut(lauf) < 150


        linescanSelectorTmp_I = linescanSelector{lauf}

       % length(linescanSelectorTmp_I);
       % midOut(lauf);

        addZero = zeros(1,150 - midOut(lauf))';

        linescanSelectorCrop = cat(1,addZero,linescanSelectorTmp_I);
        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorCrop,addZero);

        linescanSelectorCropStoreRed(:,lauf) = linescanSelectorCrop(1:300);
    elseif midOut(lauf) > 150

         linescanSelectorTmp_I = linescanSelector{lauf}


        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1,linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreRed(:,lauf) =  linescanSelectorCrop(midOut(lauf)-149:300+(midOut(lauf)-150));


    elseif midOut(lauf) == 150

        addZero =zeros(1,300)'
        linescanSelectorCrop = cat(1, linescanSelectorTmp_I,addZero);

        linescanSelectorCropStoreRed(:,lauf) = linescanSelectorCrop(1:300);

    end
end

 linescanSelectorCropStorePositionVec = [(-149:1:150)*voxelX_mum]';
linescanSelectorCropStoreTimeVec = timeVec

 FWHMOrigStoreRED =  FWHMOrigStore';
 
 FWHMTimeVecRed = timeVec
    %%%%% Fit cropped 10% above background
    %plot(x_ax,y_fitted,'r')

    %%%%% Fit non cropped curves
    %plot(x_axOrig,y_fittedOrig,'r')
     mkdir([curdir '\' 'redFit']);
       mkdir([curdir '\' 'greenFit']);

    for lauf = 20:length(linescanSelectorIn)
        redMaxTmp = maxSelector(lauf)-maxSelector(lauf);
        yGreenTmp = linescanSelector{lauf}
        h = figure(1) 
        plot(((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum),yGreenTmp,'c')
        hold on
       plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
          axis([-20 +20 0 500]) 
  %title('Gaussian fit to MRLIIC width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
   title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16);
    ylabel('Intensities [A.U.]' ,'FontSize', 16);
    legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(FWHMTimeVec(lauf)) 's' ])
        pause(0.2)
       
       
       %print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5pxAverage_5pxBF.pdf']);%tifCurvetifFilename);
   
     %   print(h,'-dpdf', [curdir '\' 'redFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5percAVG.pdf']);%tifCurvetifFilename);
   close all
    end
    

   