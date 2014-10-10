function [FWHMOrig, x_axOrig, y_fittedOrig, centeredCurve] = doFWHM_Linescan(inputLinescan,centerPosition,voxelX_mum,timeSelection) 


%centerPosition = maxSelector;
%timeSelection =    lauf
%inputLinescan=  linescanSelector{lauf}
%inputLinescan = yGreenTmp;
% centerPosition = redMax1;

   
  
      
      
      %%%%%%%%%% This part uses the entire linescan to fit the data
     % timeSelection=30;
     centeredCurve = ((0-centerPosition(timeSelection))*voxelX_mum:voxelX_mum:(length(inputLinescan)-(centerPosition(timeSelection)+1))* voxelX_mum);
    % plot(centeredCurve,inputLinescan)
     
    
% Set up fittype and options.
    zeroPos = find(centeredCurve ==0);
    %zeroValue = inputLinescan(zeroPos-10:zeroPos+10);
    xValue = inputLinescan(zeroPos-10:zeroPos+10);
    maxX = max(xValue)
    newZero = find(inputLinescan==maxX)
    newY = centeredCurve(newZero)
    
    
    ft = fittype( 'gauss1' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.Robust = 'Bisquare';
    opts.StartPoint = [maxX newY 0.5];
    opts.Upper = [Inf Inf Inf];

    success = 0;

    try
    
    % [fittedmodel, goodness, output] = fit(centeredCurve',inputLinescan,'gauss1');
      [fittedmodel, goodness, output] = fit(centeredCurve',inputLinescan,ft,opts);
     
      sigma1Orig = fittedmodel.c1/sqrt(2); %%%% Gives sigma
          FWHMOrig = 2*sqrt(2*log(2))*sigma1Orig; %%% Calculates the FWHM

        x_axOrig = -400:0.01:400; % axis
        y_fittedOrig = fittedmodel(x_axOrig); 
        success = 1;

        catch E
    end
    
    if  success == 1
    else
          FWHMOrig = 500;
           x_axOrig=0
            y_fittedOrig =0;
        
end
        %  fwhmTest1Orig = fwhm(x_axOrig,y_fittedOrig);

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

    
 