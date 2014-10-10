addpath('/Users/spira/Documents/Matlab_scripte/Image_Processing_utils')

addpath('/Users/spira/Documents/Matlab_scripte/tiffIO')
addpath('/Users/spira/Documents/Matlab_scripte/')

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;


[greenOut redOut midOut] = doAverageFlanks(redMax1,redMax2,yGreen1,yGreen2,yRed1,yRed2)

 clear linescanSelectorTmp
    linescanSelectorTmp = zeros(300,43);
     linescanSelectorTmp  = double(linescanSelectorTmp);


linescanSelectorIn = greenOut;


for lauf = 1:length(linescanSelectorIn)

   linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yGreenNonPol1{lauf})+mean(yGreenNonPol2{lauf}) )/ 2;

end

%%%%%% The following section shifts the arrays of the linescan in a way
%%%%%% that every linescan has 300px in lenght and is centered.
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

%%%%%% The following section shifts the arrays of the linescan in a way
%%%%%% that every linescan has 300px in lenght and is centered.


 clear linescanSelectorTmp
    linescanSelectorTmp = zeros(300,43);
     linescanSelectorTmp  = double(linescanSelectorTmp)
midOut


linescanSelectorIn = redOut;


for lauf = 1:length(linescanSelectorIn)

   linescanSelector{lauf} =  linescanSelectorIn{lauf}./(mean(yRedNonPol1{lauf})+mean(yRedNonPol2{lauf}) )/ 2;

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
figure(1)


imshow( linescanSelectorCropStoreGreen(:,:),[])

figure(2)

imshow( linescanSelectorCropStoreRed(:,:),[])

clear  greenStackRGKymo redStackRGBKymo imgMergeKymo

%for lauf =1 :length(linescanSelectorIn)
    greenStackRGKymo(:,:,:) = ind2rgb(normalizedImage(linescanSelectorCropStoreGreen(:,:)),green);
    redStackRGBKymo(:,:,:) = ind2rgb(normalizedImage(linescanSelectorCropStoreRed(:,:)),red);
             
    imgMergeKymo(:,:,:,:) =  greenStackRGKymo(:,:,:) +  redStackRGBKymo(:,:,:);
    figure(3)
imshow(imgMergeKymo,[])

 tiffwrite_mat(normalizedImage(linescanSelectorCropStoreGreen(:,:)), 'MyoII_Kymograph');
  tiffwrite_mat(normalizedImage(linescanSelectorCropStoreGreen(:,:)), 'SiR_Kymograph');
  
  imwrite(imgMergeKymo, 'SiR_MyoII_Kymo_double','tiff','Compression','none', ...
                                                          'resolution', [300 300]);  
                                                                        
  
  %tiffwrite_RBG(imgMergeKymo(:,:,:,:), 'SiR_Myo_Kymograph'); 
           

%for lauf = 1:length(linescanSelectorIn)
    
    %plot(1:length(linescanSelectorCropStoreGreen(:,lauf)),linescanSelectorCropStoreGreen(:,lauf))
  %  plot(1:length(linescanSelectorCropStoreRed(:,lauf)),linescanSelectorCropStoreRed(:,lauf))

  %  pause(0.2)
%end

%plot(1:length(linescanSelectorCropStoreRed(:,lauf)),linescanSelectorCropStoreRed(:,lauf))


%end