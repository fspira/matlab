%%%%% Plot cleavage furrow ingression linescan curves
FWHMStore ={};
FWHMOrigStore ={};

curdir = pwd;
%%%% Select frame to plot

%%%% Select the center position
maxSelector = redMax2;

%%%% Select the linescanchannel
linescanSelector = yGreen2;
%linescanSelector = yRed2;

lauf = 22


%%%% find frame relative to anaphase onset
time1 = timeVec(lauf);



%%%% load the data from the linescan and find the center position
redMaxTmp = maxSelector(lauf)-maxSelector(lauf);
yGreenTmp = linescanSelector{lauf}- min(linescanSelector{lauf})
%%%% Calculated the full width half maximum
[FWHM, y_fitted, x_ax, maxDistance, FWHMOrig, x_axOrig, y_fittedOrig, centeredCurve,cropedLinescan] = doFWHM_Linescan(yGreenTmp,maxSelector,voxelX_mum,lauf) ;
      
FWHMStore{1} = FWHM;
FWHMOrigStore{1} = FWHMOrig;

h=figure(1)
%%%% Plot the data
plot(((0-maxSelector(lauf))*voxelX_mum:voxelX_mum:(length(yGreenTmp)-(maxSelector(lauf)+1))* voxelX_mum),yGreenTmp,'c')

%legend([num2str(time1) ' s'])
hold on

%%%%% Fit cropped 10% above background
%plot(x_ax,y_fitted,'r')

%%%%% Fit non cropped curves
%plot(x_axOrig,y_fittedOrig,'r')


%%%% Select frame to plot

lauf =27
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

lauf = 31
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
axis([-10 +10 0 800]) 
    yL = get(gca,'YLim');
    line([redMaxTmp redMaxTmp],yL,'Color','r'); 
    line([redMaxTmp-(30*voxelX_mum) redMaxTmp-(30*voxelX_mum)],yL,'Color','g');
     line([redMaxTmp+(30*voxelX_mum) redMaxTmp+(30*voxelX_mum)],yL,'Color','g');
   xlabel ('Distance [µm]','FontSize', 20);
   ylabel('Intensity [A.U.]','FontSize', 20);
   title(['SiR-actin RAW Intensity Sliding Window:' tifFilename],'FontSize', 20);
   

   print(h,'-dpdf', [curdir '/' tifFilename,'_FWHM_CleavageFurrow_MLCIIB.pdf']);%tifCurvetifFilename);


   