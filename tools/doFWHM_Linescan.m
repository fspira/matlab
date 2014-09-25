function [FWHM, y_fitted, x_ax, maxDistance, FWHMOrig, x_axOrig, y_fittedOrig, centeredCurve,LinescanAnalysis] = doFWHM_Linescan(inputLinescan,centerPosition,voxelX_mum,timeSelection) 


   lauf = 1;

%inputLinescan = yGreenTmp;
% centerPosition = redMax1;

    LinescanTmp =  inputLinescan;
    LinescanBkg = (mean(LinescanTmp(1:10))+mean(LinescanTmp(end-10:end)))/2;
       centerPixel =round(length(LinescanTmp)/2);
    %%%%%% Add 10% to Bkg to ensure that the signal chosen for analysis is
    %%%%%% real
    LinescanBkgSub = LinescanTmp - (LinescanBkg+((LinescanBkg/100)*20));
    
    %plot(1:length(LinescanBkgSub), LinescanBkgSub)
    %hold on
   % plot(1:length(inputLinescan), inputLinescan)
    
    %%%%%% Eliminate elements smaller than zero
    LinescanBkgSubTmp = LinescanBkgSub >= 0;
   
    runMarker = 1;
    clear leftEdge;
    clear rightEdge;
    
    for subrun = 1:centerPixel-1
        
      if  LinescanBkgSubTmp(centerPixel-subrun) == 1
        
      else
            leftEdge(runMarker) = centerPixel-subrun
            runMarker = runMarker +1;
      end
            
    end
      runMarker = 1;
     for subrun = 1:  (length(LinescanBkgSubTmp) - centerPixel)-1
        
       if LinescanBkgSubTmp(centerPixel+subrun) == 1
        
       else
            rightEdge(runMarker) = centerPixel+subrun;
             runMarker = runMarker +1;
       end
            
     end
    
     LinescanAnalysis = LinescanBkgSub( leftEdge(1):rightEdge(1));
     
    
     
     %plot(1:length(LinescanAnalysis), LinescanAnalysis);
    
    
      
      %%%%% Calculate euclidian distance - sum the distance between the
      %%%%% center of each pixel
      
        furrowAnalysis = length(LinescanBkgSub( leftEdge(1):rightEdge(1)));
             clear distTmp
            for lauf_1 = 1:furrowAnalysis-1
                distTmp(lauf_1) = pdist2(LinescanAnalysis(lauf_1,:),LinescanAnalysis(lauf_1+1,:));


            end
     
            
             euclidianDistanceFurrow(lauf) = sum(distTmp);
      
      
           distance = 0+voxelX_mum:voxelX_mum: length(LinescanAnalysis) * voxelX_mum;
           distanceOriginal = 0+voxelX_mum:voxelX_mum: length(inputLinescan) * voxelX_mum;
           maxDistance = length(LinescanAnalysis) * voxelX_mum;
           maxDistanceOriginal = length(inputLinescan) * voxelX_mum;

         [fittedmodel, goodness, output] = fit(distance',LinescanAnalysis,'gauss1');
        % [fittedmodel, goodness, output] = fit(distanceOriginal',inputLinescan,'gauss1');
         
          sigma1 = fittedmodel.c1/sqrt(2); %%%% Gives sigma
          FWHM = 2*sqrt(2*log(2))*sigma1; %%% Calculates the FWHM

        x_ax = -400:0.01:400; % axis
        y_fitted = fittedmodel(x_ax); 


          fwhmTest1 = fwhm(x_ax,y_fitted);

      %  h = figure

       % plot(distance,LinescanAnalysis,'r');
     %  plot(distanceOriginal, inputLinescan,'b')
      %  hold on
      %  plot(x_ax,y_fitted)
      %  axis([-4 14 0 max(LinescanAnalysis)+40])  ;


       % title(['FWHM', char(lauf)])
       % xlabel ('Distance [µm]','FontSize', 16);
       %  ylabel('Intensities [A.U.]','FontSize', 16);



       % print(h,'-dpdf', [curdir '/' tifFilename,'_FWHM.pdf']);%tifCurvetifFilename);
     %   close all

         
      FWHM_Store{lauf} = FWHM;
      y_fitted_Store{lauf} = y_fitted;
      x_ax_Store{lauf} = x_ax;
      maxDistance_Store{lauf} = maxDistance;
      
  
      
      
      %%%%%%%%%% This part uses the entire linescan to fit the data
     % timeSelection=30;
     centeredCurve = ((0-centerPosition(timeSelection))*voxelX_mum:voxelX_mum:(length(inputLinescan)-(centerPosition(timeSelection)+1))* voxelX_mum);
    % plot(centeredCurve,inputLinescan)
     
     [fittedmodel, goodness, output] = fit(centeredCurve',inputLinescan,'gauss1');
     
      sigma1Orig = fittedmodel.c1/sqrt(2); %%%% Gives sigma
          FWHMOrig = 2*sqrt(2*log(2))*sigma1Orig; %%% Calculates the FWHM

        x_axOrig = -400:0.01:400; % axis
        y_fittedOrig = fittedmodel(x_axOrig); 


          fwhmTest1Orig = fwhm(x_axOrig,y_fittedOrig);

      %  h = figure

       % plot(distance,LinescanAnalysis,'r');
      % plot(centeredCurve, inputLinescan,'b')
       % hold on
       % plot(x_axOrig,y_fittedOrig)
       % axis([-14 14 0 max(inputLinescan)+40])  


       % title(['FWHM input linescan', char(lauf)])
       % xlabel ('Distance [µm]','FontSize', 16);
       % ylabel('Intensities [A.U.]','FontSize', 16);
     
     
      
 % provides all the fitted values for  
                        %%%%the specified x-axis

    
 