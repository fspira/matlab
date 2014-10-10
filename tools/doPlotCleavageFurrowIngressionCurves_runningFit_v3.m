%%%%% Plot cleavage furrow ingression linescan curves


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
FWHMOrigStore =[];

%%%% Average flanks
[greenOut redOut midOut] = doAverageFlanks(redMax1,redMax2,yGreen1,yGreen2,yRed1,yRed2)


curdir = pwd;
%%%% Select frame to plot

%%%% Select the center position
maxSelector = midOut;

%%%% Select the linescanchannel
%linescanSelectorIn = greenOut;
%linescanSelectorIn = redOut;

%%%% Use the poles as a reference for normalization
%%%% Use the total intensity for normalization

%for lauf = 1:length(linescanSelectorIn)

 %   linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yGreenNonPol1{lauf})+mean(yGreenNonPol2{lauf}) + mean(yGreen1{lauf}) + mean(yGreen2{lauf}))/4;

%end


%%% Smooth the curve a sliding average of 5 is used
%%% Since the window creates artifacts at the left edge, a dummy value is
%%% first attached to the left edge and (mean of 5 Values and subsequently
%%% removed

windowSize = 5;

b = (1/windowSize)*ones(1,windowSize)
a = 1;

for lauf = 1:length(linescanSelectorIn)
    linescanTmp = linescanSelectorIn{lauf};
    linescanTmpMean = mean(linescanTmp(1:5))
    linescanFill(1:5) = linescanTmpMean ;
    linescanTmp = cat(1,linescanFill',linescanTmp);
    
    linescanSelectorTmp = filter(b,a,linescanTmp);
    linescanSelectorIn{lauf} = linescanSelectorTmp(6:length(linescanSelectorTmp))
    % plot(length(linescanSelectorIn{lauf}),linescanSelectorIn{lauf},'b')
end

%%%%% use a flatline background correction - the left and right 20 pixels
%%%%% are used for correction

for lauf = 1:length(linescanSelectorIn)
     linescanTmp = linescanSelectorIn{lauf};
    [y,yfit] = bf( linescanTmp,[1:5,length(linescanTmp)-5:length(linescanTmp)],'linear');
    
    y_Indx = y > 0;
    yNew = y_Indx .* y;

    linescanSelectorIn{lauf} = yNew;
    
end

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

for lauf = 1:length(linescanSelectorIn)

   linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yGreenNonPol1{lauf})+mean(yGreenNonPol2{lauf}) )/ 2;
   
   %linescanSelector{lauf} = linescanSelector{lauf} - min(linescanSelector{lauf}); 
   
   
end


for lauf = 1:length(linescanSelectorIn)

   linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yRedNonPol1{lauf})+mean(yRedNonPol2{lauf}) )/ 2;
  % linescanSelector{lauf} = linescanSelector{lauf} - min(linescanSelector{lauf}); 
   
   
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

 FWHMOrigStore =  FWHMOrigStore';
 
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
          axis([-20 +20 0 1]) 
  title('Gaussian fit to MyoII width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16);
    ylabel('Intensities [A.U.]' ,'FontSize', 16);
    legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(FWHMTimeVec(lauf)) 's' ])
        pause(0.2)
       
       
       print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5pxAverage_5pxBF.pdf']);%tifCurvetifFilename);
   
       % print(h,'-dpdf', [curdir '\' 'redFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_noAverage_NopxBFminAvgSub.pdf']);%tifCurvetifFilename);
   close all
    end
    
    plot(1:length(FWHMOrigStore),FWHMOrigStore,'.')
    axis([15 40 0 25])
    %%%% Select frame to plot

    lauf =32
    %%%% find frame relative to anaphase onset
    time2 = timeVec(lauf);

    %%%% load the data from the linescan and find the center position
    redMaxTmp = maxSelector(lauf)-maxSelector(lauf);
    yGreenTmp = linescanSelector{lauf}- min(linescanSelector{lauf})



    %%%% Calculated the full width half maximum
    [FWHM, y_fitted, x_ax, maxDistance, FWHMOrig, x_axOrig, y_fittedOrig, centeredCurve,cropedLinescan] = doFWHM_Linescan(yGreenTmp,maxSelector,voxelX_mum,lauf) ;

    FWHMStore{2} = FWHM;
    FWHMOrigStore{2} = FWHMOrig;

    %%%% Plot the data
    plot(((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum),yGreenTmp,'b')
    %legend([num2str(time2) ' s'])
    %%%%% Fit non cropped curves
    %plot(x_axOrig,y_fittedOrig,'r')


    %%%%% Fit cropped 10% above background
    %plot(x_ax,y_fitted,'r')



    %%%% Select frame to plot

    lauf = 36
    %%%% find frame relative to anaphase onset
    time3 = timeVec(lauf);
    %%%% load the data from the linescan and find the center position
    redMaxTmp = maxSelector(lauf)-maxSelector(lauf);
    yGreenTmp = linescanSelector{lauf}- min(linescanSelector{lauf})

    %%%% Calculated the full width half maximum
    [FWHM, y_fitted, x_ax, maxDistance, FWHMOrig, x_axOrig, y_fittedOrig, centeredCurve,cropedLinescan] = doFWHM_Linescan(yGreenTmp,maxSelector,voxelX_mum,lauf) ;

    FWHMStore{3} = FWHM;
    FWHMOrigStore{3} = FWHMOrig;


    %%%% Plot the data
    plot(((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum),yGreenTmp,'r')

    %%%%% Fit cropped 10% above background
    %plot(x_ax,y_fitted,'r')

    %%%%% Fit non cropped curves
    %plot(x_axOrig,y_fittedOrig,'r')


    legend([num2str(time1) ' s'],[num2str(time2) ' s'],[num2str(time3) ' s'])
    %legend([num2str(time1) 's'],['FWHM ' num2str(FWHMStore{1}) ' s'],[num2str(time2) 's'],['FWHM '  num2str(FWHMStore{2}) ' s'],[num2str(time3) ' s'],['FWHM ' num2str( num2str(FWHMStore{3})) 's'])

    %%%% label axis
    axis([-10 +10 0 1]) 
        yL = get(gca,'YLim');
        line([redMaxTmp redMaxTmp],yL,'Color','r'); 
        line([redMaxTmp-(30*voxelX_mum) redMaxTmp-(30*voxelX_mum)],yL,'Color','g');
         line([redMaxTmp+(30*voxelX_mum) redMaxTmp+(30*voxelX_mum)],yL,'Color','g');
       xlabel ('Distance [µm]','FontSize', 20);
       ylabel('Intensity [A.U.]','FontSize', 20);
       title(['SiR-actin RAW Intensity Sliding Window:' tifFilename],'FontSize', 20);
end

   print(h,'-dpdf', [curdir '/' tifFilename,'_FWHM_CleavageFurrow_MLCIIB.pdf']);%tifCurvetifFilename);


   