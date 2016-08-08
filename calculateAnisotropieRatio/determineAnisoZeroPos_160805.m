
%%%%%
%%%%% Anisotropy correction model
%%%%%

load('/Users/spira/Desktop/3720/160713/metaphaseNormalization/stressFibers_160726.mat')
load('/Users/spira/Desktop/3720/160713/metaphaseNormalization/fitresultMetaphase_160726.mat')
load('/Users/spira/Desktop/3720/160713/metaphaseNormalization/fitresultStressFibers_160726.mat')


%load('/Users/spira/Desktop/3720/cellPoolForCorrection/metaphaseNormalization/fitresultStressFibers_3720_160725.mat')

counter = -35
for lauf = 1:360
    
    
   fitValues(lauf) = feval(fitresultstressFibers,lauf);
   fitValuesMeta(lauf) = feval(fitresultMetaphase,counter);
    counter = counter+1
    
end
h = figure(1)

plot(fitValues,'b')
hold on
plot(fitValuesMeta,'r')
axis([0 180 0.7 2.2])
plot([0 180], [meanStress meanStress],'--b')
  plot([0 180], [meanMeta meanMeta],'--r')
  
xlabel ('Angle (deg)','FontSize', 16);
ylabel('Fluorescence anisotropy','FontSize', 16);
title(['Raw ratio'],'FontSize', 16);


legend('Anisotropy stress fibers', ...
  'Anisotropy metaphase',...
  'Mean stress fibers',...
  'Mean metaphase')
print(h,'-dpdf', ['rawRatio.pdf']);%tifCurvetifFilename);


meanStress = mean(fitValues)
meanMeta = mean(fitValuesMeta)
metaOffset = meanStress - meanMeta
%%%%%%%
%%%%%%% mean normalization
%%%%%%%
fitValuesMetaMeanNorm = fitValuesMeta ./ meanMeta
fitValuesMeanNorm = fitValues ./ meanStress


maxStressMeanNorm = max(fitValuesMeanNorm)
maxMetaMeanNorm = max(fitValuesMetaMeanNorm)

normMeta = (fitValuesMetaMeanNorm - minNormMeta) ./ (maxStressMeanNorm -minNormMeta)
max(normMeta)
minNormMeta = min(fitValuesMetaMeanNorm)
normMetaCorr = normMeta - minNormMeta

max(normMetaCorr)

h = figure(2)
plot(normMeta)
hold on
plot(normMetaCorr)


plot(fitValuesMeanNorm,'b')
hold on
plot(fitValuesMetaMeanNorm,'r')
plot([0 180], [1 1],'k')
axis([0 180 0.5 1.5])
  
    
xlabel ('Angle (deg)','FontSize', 16);
ylabel('Fluorescence anisotropy','FontSize', 16);
title(['Mean normalized anisotropy'],'FontSize', 16);

legend('Anisotropy stress fibers', ...
  'Anisotropy metaphase', ...
  'Isotropic boundary')

print(h,'-dpdf', ['_meanNormalizedanisotropy.pdf']);%tifCurvetifFilename);

fitValuesMetaOffsetCorr = fitValuesMeta + metaOffset

plot(fitValues)
hold on
plot(fitValuesMetaOffsetCorr)

maxMetaCorr = max(fitValuesMetaOffsetCorr)
maxStress = max(fitValues)
minStress = min(fitValues)
anisoStress = mean(fitValues)

normMetaMaxAniso = (maxMetaCorr - anisoStress) ./ (maxStress - anisoStress)
normMetaMaxMin = (maxMetaCorr -minStress) ./ (maxStress - minStress)

for lauf =1:180
    
fitMetaValuesNorm(lauf) = fitValuesMetaMeanNorm(lauf) ./ fitValuesMeanNorm(lauf)
fitMetaValuesMetaNorm(lauf) = fitValuesMetaMeanNorm(lauf) ./ fitValuesMetaMeanNorm(lauf)

end

maxStressNorm = max(fitValuesMeanNorm)
minStressNorm = min(fitValuesMeanNorm)
anisoStressNorm = mean(fitValuesMeanNorm)



maxMetaNorm = max(fitValuesMetaMeanNorm)
minMetaNorm = min(fitValuesMetaMeanNorm)
%anisoStressNorm = mean(fitMetaValuesMetaNorm)

figure(2)
plot(fitValuesMeanNorm)
hold on
plot(fitValuesMetaMeanNorm)

h = figure(2)


plot(fitValuesMeanNorm,'b')
hold on
plot(fitValuesMetaMeanNorm,'r')
%plot([0 180], [1 1],'k')
axis([0 180 0.5 1.5])
  
plot(fitMetaValuesNorm,'--b')
hold on
plot(fitMetaValuesMetaNorm,'--r')

    
xlabel ('Angle (deg)','FontSize', 16);
ylabel('Fluorescence anisotropy','FontSize', 16);
title(['Mean normalized anisotropy'],'FontSize', 16);

legend('Anisotropy stress fibers', ...
  'Anisotropy metaphase', ...
  'Stress model for correction', ...
  'Metaphase model for corrected')

print(h,'-dpdf', ['CorrectedNormalizedAnisotropyI.pdf']);%tifCurvetifFilename);



axis([0 180 0.5 1.5])
%%%%% normalize meta values 


maxStressNorm = max(fitValuesMeanNorm)
minStressNorm = min(fitValuesMeanNorm)
anisoStressNorm = mean(fitValuesMeanNorm)



maxMetaNorm = max(fitValuesMetaMeanNorm)
minMetaNorm = min(fitValuesMetaMeanNorm)

metaNormStress_I = (fitMetaValuesMetaNorm - minStressNorm) ./ (maxStressNorm - minStressNorm)

metaNormStress_II = (maxMetaNorm - 1) ./ (maxStressNorm - 1)

metaNormStress_III = (maxMetaNorm - minMetaNorm) ./ (maxStressNorm - minMetaNorm)

metaNormStress_IV= (maxMetaNorm - minStressNorm) ./ (maxStressNorm - minStressNorm)

metaNormStress_V = (fitMetaValuesMetaNorm - minMetaNorm) ./ (maxStressNorm - minMetaNorm)

%metaNormStress_VI = (fitMetaValuesMetaNorm - minStressNorm) ./ (maxStressNorm - minStressNorm)

h = figure(3)

plot(metaNormStress_I,'k')
hold on

plot([0 180], [metaNormStress_II,metaNormStress_II],'b')

plot([0 180],[metaNormStress_III,metaNormStress_III],'r')

plot([0 180],[metaNormStress_IV,metaNormStress_IV],'g')
plot(metaNormStress_V,'c')



axis([0 180 0 1])



legend('Norm values between min and max stress fibers', ...
  'Max meta between aniso and max stress', ...
  'Max meta between min meta and max stress', ...
  'Max meta between min stress and max stress', ...
  'Norm values between min meta and max stress')



xlabel ('Angle (deg)','FontSize', 16);
ylabel('Normalized anisotropy','FontSize', 16);
title(['Normalization methods'],'FontSize', 16);


print(h,'-dpdf', ['Normalization methods.pdf']);%tifCurvetifFilename);


metaNormMeta = (fitValuesMetaMeanNorm - 1) ./ (maxStressNorm - 1)
max(metaNormMeta)



%%%%%%
%%%%%% normalize meta values by moving the 
%%%%%%

plot(metaNorm)



    fitValuesMetaNorm = fitValues  ./ mean(fitValues);
       
    [fitresultstressFibersNorm, gofStress] = createFitSineSquareFixPhase([1:180], fitValuesMetaNorm')
    
         
    
xlabel ('Time [s]','FontSize', 16);
ylabel('Ratio [A.U.]','FontSize', 16);
title(['Ratio green/red furrow and pole' tifFilename],'FontSize', 16);

print(h,'-dpdf', [tifFilename,'_ratioFurrowPole.pdf']);%tifCurvetifFilename);


