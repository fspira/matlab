%%%%%% Read and plot merged Data from furrow Ingression experiments

clear all

%diameter1 = xlsread('mergeDistances','Sheet1');
[~, ~, raw] = xlsread('/Users/spira/Desktop/Rpe1_siRActin_H2BRFP/analysis/mergeDistancesAutoH2B.xlsx');
raw = raw(4:end,:);

%%%%%% Create output variable
diameter1 = cell2mat(raw);
clear raw



%diameter1 = xlsread('mergeDistances','Sheet1');
[~, ~, raw] = xlsread('/Users/spira/Desktop/Rpe1_siRActin_H2BRFP/analysis/mergeDistancesTest.xlsx');
raw = raw(4:end,:);

%%%%%% Create output variable
testVariables = cell2mat(raw);
clear raw

[m n p] = size(diameter1);

time = diameter1(:,21);
timeFit = diameter1(:,1:5:end-10)
contractileRingDiameter = diameter1(:,2:5:end-10)
polePoleDistance = diameter1(:,3:5:end-10)
chromatinChromatinDistance = diameter1(:,4:5:end-10)
chromatinPoleDistance = diameter1(:,5:5:end-10)

%timeTest = diameter1(:,21);
timeFitTest = diameter1(:,end-9:5:end)
contractileRingDiameterTest = diameter1(:,end-8:5:end)
polePoleDistanceTest = diameter1(:,end-7:5:end)
chromatinChromatinDistanceTest = diameter1(:,end-6:5:end)
chromatinPoleDistanceTest = diameter1(:,end-5:5:end)


timeTest = testVariables(:,1)
ringDiameterTest = testVariables(:,2)
chromatinDistanceTest=  testVariables(:,4)



%%%% Cell number 5 is bigger (double the size, probably polyploid and will
%%%% be removed. Chromatin at position 4 fails to correctly divide
%%%% neglected for training set

[mmm nnn ppp] = size(chromatinChromatinDistance)

timeFitTmp = timeFit(:,1:3);timeFit= cat(2,timeFitTmp,timeFit(:,6:nnn));
timeFit = cat(2,timeFit,timeFitTest(:,1))


contractileRingDiameterTmp = contractileRingDiameter(:,1:3);
contractileRingDiameter = cat(2,contractileRingDiameterTmp,contractileRingDiameter(:,6:nnn))
contractileRingDiameter = cat(2,contractileRingDiameter, contractileRingDiameterTest(:,1))

chromatinChromatinDistanceTmp = chromatinChromatinDistance(:,1:3)
chromatinChromatinDistance = cat(2,chromatinChromatinDistanceTmp,chromatinChromatinDistance(:,6:nnn))
chromatinChromatinDistance = cat(2,chromatinChromatinDistance,chromatinChromatinDistanceTest(:,1))

[mmm nnn ppp] = size(chromatinChromatinDistance)
 timeFitCell = {};
 chromatinChromatinDistanceCell = {}
 contractileRingDiameterCell = {}

for lauf = 1:nnn
    chromatinChromatinDistanceTmp = chromatinChromatinDistance(:,lauf)
    chromatinChromatinDistanceTmp1 = find(chromatinChromatinDistanceTmp > 0)
    matrixStart = chromatinChromatinDistanceTmp1(1);
    matrixEnd =  chromatinChromatinDistanceTmp1(end);
    
    timeFitCell{lauf}= timeFit(matrixStart:matrixEnd,lauf)
    chromatinChromatinDistanceCell{lauf} = chromatinChromatinDistance(matrixStart:matrixEnd,lauf);
    contractileRingDiameterCell{lauf} = contractileRingDiameter(matrixStart:matrixEnd,lauf)
end


   h= figure; 
   
  subplot(2,1,1) 
   hold on
for lauf = 1:nnn
   plot(timeFitCell{:,lauf},contractileRingDiameterCell{:,lauf})
   pause(0.2)
end
    hold on
    axis([0 500 0 20])  
    title('H2BCherry sirActin - contractile ring diameter','FontSize', 16)
    xlabel ('Time [s]','FontSize', 16);
    ylabel('Distance [\mum]' ,'FontSize', 16);

 
           
    %%%%% plot chromatin-chromatin distance
    
    
    subplot(2,1,2)
    hold on
   for lauf = 1:nnn
    plot(timeFitCell{:,lauf},chromatinChromatinDistanceCell{:,lauf})
    pause(0.2)
    
   end
    hold on
    axis([0 500 0 20])  
    title('H2BCherry sirActin - chromatin chromatin distance','FontSize', 16)
    xlabel ('Time [s]','FontSize', 16);
    ylabel('Distance [\mum]' ,'FontSize', 16);
    
    

    
    %%%% Normalize data to anaphase Onset - as determined by chromatin
    %%%% segregation
    
    
    
 orient landscape;
     print(h,'-dpdf', ['DistancesSingleTracksAnalysis.pdf']);%tifCurvetifFilename);
    close all
           
    
    
  %  [m n p] = size(greenStore);

[stdCRDiameter,meanCRDiameter] = calculateStd(contractileRingDiameter(:,:))
[stdPolePole,meanPolePole] = calculateStd(polePoleDistance)
[stdChrChr,meanChrChr] = calculateStd(chromatinChromatinDistance(:,:))
[stdChrPole,meanChrPole] = calculateStd(chromatinPoleDistance)





h1 = errorbar(time,meanCRDiameter,stdCRDiameter,'r')
hold on    
h2 = errorbar(time,meanChrChr,stdChrChr,'b')

 title('Chormatin distance and contractile ring diameter','FontSize', 16)
    xlabel ('Time [s]','FontSize', 16);
    ylabel('Distance [µm]' ,'FontSize', 16);

    
%legend([h1 h2], ...
%  'Contractile ring diameter', ...
%   'Chromatin-chromatin distance')


%%%%%%%%%%%%% plot cortical intensities

 

h=  figure(2)
   
    hold on
   
    %errorbar(time,greenMean,stdGreenLine,'g')
   errorbar(time,meanFurrowNorm,stdFurrowNorm,'r')
   
    
    title('H2Bcherry sirActin - pole intensity norm to ana','FontSize', 20)
    xlabel ('Time [s]','FontSize', 16);
    ylabel('Rel. Intensity [A.U.]' ,'FontSize', 16);


 orient landscape;
     print(h,'-dpdf', ['ContractileRingDiameter.pdf']);%tifCurvetifFilename);
    close all
    [stdCRDiameter,meanCRDiameter] = calculateStd(contractileRingDiameter)
[stdPolePole,meanPolePole] = calculateStd(polePoleDistance)
[stdChrChr,meanChrChr] = calculateStd(chromatinChromatinDistance)
[stdChrPole,meanChrPole] = calculateStd(chromatinPoleDistance)

    
x = time;
y1 = meanFurrowNorm;
y2 = meanCRDiameter;
y3 = meanPolePole;
y4 =meanChrChr;
y5 =meanChrPole;


[AX,H1,H2] = plotyy(x,y1,x,y2, 'plot');
  
  set(get(AX(2),'Ylabel'),'String','Distance [µm]','FontSize', 16) 
 %  set(get(AX(1),'Ylabel'),'String','Distance [µm]','FontSize', 16) 
  set(get(AX(1),'Ylabel'),'String','Normalized Mean Intensity [A.U.]','FontSize', 16) 
       
hold(AX(1), 'on')
hold(AX(2), 'on')

% Plot the third curve
h3 = plot(x, y4, 'c', 'Parent', AX(2));
    

% Plot the third curve
%h4 = plot(x, y5, 'r', 'Parent', AX(2));
%h3 = errorbar(time,meanAnillinH2B_Chrom,stdAnillinH2B_Chrom,'g')

%h4 = plot(x, y3, 'r', 'Parent', AX(2));
%h4 = errorbar(time,meanAnillinH2B_Ring,stdAnillinH2B_Ring,'g')


%h3 = errorbar(time,meanLifeactH2B_Chrom,stdLifeactH2B_Chrom,'g')

legend([H1 H2 h3], ...
  'Mean furrow intensity', ...
  'Contractile ring diameter', ...
  'Chromosome-chromosome Distance')%, ...
 % 'Chromosome-pole Distance')


 title('sirActin-H2BCherry - chromosome cortex distances and intensity','FontSize', 16)
    xlabel ('Time [s]','FontSize', 16);
  
 %%%%%%% This part generates a scatterplot between contractile ring
 %%%%%%% diameter and H2B distance.
 
 
 h=figure(1)
 hold on
 for lauf =1:6
 plot( chromatinChromatinDistanceCell{:,lauf},contractileRingDiameterCell{:,lauf},'x')
 
% plot3(chromatinChromatinDistanceCell{:,lauf},contractileRingDiameterCell{:,lauf},timeFitCell{:,lauf},'x')
 end
 h= figure(1) 
 hold on
 for lauf = 1:6
     timeDummy = length(timeFit(:,lauf))
    scatter((chromatinChromatinDistanceCell{:,lauf}), (contractileRingDiameterCell{:,lauf}),timeFitCell{:,lauf})
    
 end
 
 

 title('RPE1 H2B - siRActin - Staging. Ring dimater encodes time information - small to big','FontSize', 16)
    xlabel ('Chromatin-chromatin distance [µm]','FontSize', 16);
    ylabel ('Contractile ring diameter [µm]','FontSize', 16);
  
  

     print(h,'-dpdf', ['Scatter_ChromatinChromatin_Ring_Timeauto.pdf']);%tifCurvetifFilename);
    close all
    
    
    
  

  
% timeMerge = cat(1,time(5:end-3,1),time(1:end-3,1))
% timeMerge = cat(1,timeMerge,time(1:end-3,1))
% timeMerge = cat(1,timeMerge,time(1:end-5,1))
 
 
x1 = meanCRDiameter(1:end-3)';
x2 = meanChrChr(1:end-3)';



timeMean = time(5:end-3);
%%%% Building the data input vectors to train the model
xx1 = chromatinChromatinMerge;
xx2 =  contractileDistMerge;
x_KDE = [xx1 xx2]
x1 = []
x2 = []
t = []
for lauf = 1:6
    x1 = cat(1,x1, chromatinChromatinDistanceCell{:,lauf});
    x2 = cat(1,x2, contractileRingDiameterCell{:,lauf});
    t = cat(1,t, timeFitCell{:,lauf});
end
    

contractileRingDiameterTest(:,1)


for subrun = 1:2
    
    timeFitTestTmp = timeFitTest(11:end,subrun);
    chromatinChromatinDistanceTestTmp = chromatinChromatinDistanceTest(11:end,subrun);
    contractileRingDiameterTestSingleTmp = contractileRingDiameterTest(11:end,subrun);
    
    chromatinChromatinDistanceTestSingle{subrun} = chromatinChromatinDistanceTestTmp(isfinite(chromatinChromatinDistanceTestTmp(:, 1)), :)
     
    contractileRingDiameterTestSingle{subrun} = contractileRingDiameterTestSingleTmp(isfinite(contractileRingDiameterTestSingleTmp(:, 1)), :)
    timeFitTestSingle{subrun} = timeFitTestTmp(isfinite(timeFitTestTmp(:,1)),:)
      
    
end
 
%%%%% Organizing the test data into a variable

  %  timeToPredict = [chromatinDistanceTest(4:end,1) ringDiameterTest(4:end,1) timeTest(4:end,1)]
  timeToPredict = [  chromatinChromatinDistanceTestSingle{2}  contractileRingDiameterTestSingle{2} timeFitTestSingle{2}]
  [mm nn pp] = size(timeToPredict)
     
     
for k=4:4 
    h= figure; 
    
    %%%%% Building the model by fitting a polynomial function to the
    %%%%% trainig data
   
    m = polyfitn([x1, x2], t, k);
    
    %%%%% Test the model by evaluating test data. The test data is real
    %%%%% measured data
    %%%%% Predicted time line is the fitted curve. Data points for staging
    %%%%% will be compared against this time line.
    predicted_timeLine = polyvaln(m, [x1, x2]); 
  
    
    %%%%% predicted_time is the predicted time on basis of the learned
    %%%%% model
    predicted_time = polyvaln(m,timeToPredict(:,1:2))
             
     testArray = cat(2,predicted_time,timeToPredict(:,3))
   
    plot3(timeToPredict(:,1),timeToPredict(:,2), timeToPredict(:,3), 'ro'); hold on;
     
    plot3(timeToPredict(:,1), timeToPredict(:,2), predicted_time, 'bx')
    
    
    residuals{k} = timeToPredict(:,3) - predicted_time;
    SSresid(k) = sum(residuals{k}.^2);
    
   
    
    plot(1:k, SSresid(1:k),'x');hold on;
    %axis([0 21 0 1000000])
    
    
  
    xlabel('Contractile ring diameter [µm]')
    ylabel('Chromatin-chromatin distance [µm]')
    
    title([num2str(k),' degree']);
    hold off
  % print(h,'-dpdf', [num2str(k),' degree_TrainingModelTestData.pdf']);%tifCurvetifFilename);
  %  close all
    
         h = figure
    plot3(x1, x2, t, 'r.'); hold on; 
    plot3(x1,x2, predicted_timeLine, 'bo');
   % plot3(timeToPredict(:,1), timeToPredict(:,2), predicted_time, 'bx')
    
    xlabel('Contractile ring diameter [µm]')
    ylabel('Chromatin-chromatin distance [µm]')
    
    title([num2str(k),' degree']);
    hold off
   % print(h,'-dpdf', [num2str(k),' degree_TrainingModelRawDataTestData.pdf']);%tifCurvetifFilename);
   % close all
    
    
    h= figure
    
    plot(1:mm,residuals{k},'xb'); hold on;
    
    
    xlabel('Data point')
    ylabel('Time difference [s]')
    
    title([num2str(k),' degree - residuals']);
    hold off

   % print(h,'-dpdf', [num2str(k),' degree_Residuals.pdf']);%tifCurvetifFilename);
    %close all
    
    
    
    
    
end

%%%%% best fit with k= 4
timeToPredict = uint8([12 5])

[density flagOut predicted_time] = doMapTime(timeToPredict);

predicted_time;

flagOut



%k=9

[bandwidth,density,X,Y]=kde2d([x1 x2], 256, [0 0], [30 30]);

h = figure
%contour3(X,Y,density,50), view([0,90]),hold on
contour3(X,Y,density,50), hold on

for lauf = 1:5
plot(chromatinChromatinDistanceCell{:,lauf},contractileRingDiameterCell{:,lauf},'r.','MarkerSize',5)
end
plot(timeToPredict(:,1),timeToPredict(:,2),'x','MarkerSize',5)

    ylabel('Contractile ring diameter [µm]')
    xlabel('Chromatin-chromatin distance [µm]')
    title('Kernel density function raw data');
     print(h,'-dpdf', ['Kernel_density_function_raw_data.pdf']);
    

%[xx, yy] = meshgrid(0:0.1:20, 0:0.1:20);
%res = polyvaln(m, [xx(:),yy(:)]);
%a = reshape(predicted_time, 201, 201);
%imshow(a,[])

%%%%% The following small functions map the measured data into the density space.
%%%%% The kernel density space ranges from 1 to 256, while the X Y Values
%%%%% range from 0 to 30. The Threshold was manually determined by
%%%%% using the high value peack of the histogram as a threshold. Density
%%%%% is the computed probability space and any new point can be compared
%%%%% to this space.


ff = @(xxx) xxx/30*255+1;
density(ff(timeToPredict(:,1)),ff(timeToPredict(:,2))) > 0.004;

df = @(xx1,yy1) density(ff(xx1),ff(yy1)) > 0.004;
df(timeToPredict(1),timeToPredict(2))

 