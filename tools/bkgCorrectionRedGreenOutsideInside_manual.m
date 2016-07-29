function [BackgroundA,BackgroundB,Background1, Background2, Background3, Background4,BWHor_Out] = bkgCorrectionRedGreenOutsideInside(img1,img2,voxelX_mum)

%img1 = greenStack;

imgNorm = normalizedImage(img2(:,:,1));

[mf nf pf] = size(img2);

fh = figure(1);
 title('Mark the Background')
frapBkgIn = roipoly(imgNorm(:,:,1));
close(fh);


for lauf = 1:pf
    
    imgNorm1 = normalizedImage(img1(:,:,lauf));

imshow(imgNorm1(:,:,1),'InitialMagnification',1000)


%colormap(jet)
[xx,yy] = ginput(1)
close all

    for subrun = 1:1

    %%%% This command draws an hourglass
    %yTest = [y-10 y-10 y+10 y+10]
    %xTest = [x+10 x-10 x+10 x-10];

    %%%%%%%% setup of the ROI - 1.5µm x 4µm

    x = xx(subrun);
    y = yy(subrun);

    shortAxis = round((2 /voxelX_mum)/2);
    longAxis = round((2 / voxelX_mum)/2);


    %%%%%% Edges of the horizontal box
    yTest =    [y-shortAxis y+shortAxis y+shortAxis y-shortAxis y-shortAxis]'
    xTest =    [x-longAxis x-longAxis x+longAxis x+longAxis x-longAxis]'

    yTestVer = [y-longAxis y-longAxis y+longAxis y+longAxis y-longAxis]'
    xTestVer = [x-shortAxis x+shortAxis x+shortAxis x-shortAxis x-shortAxis]'


    %%%%%% Total number of pixels within the box
    BW_Box = (yTest(2)-yTest(1)) * (xTest(3) - xTest(2));

    %%%%%% Draw the horizontal box
    BWHor = roipoly(imgNorm1,xTest, yTest);

    %%%%% Draw the vertical box
    %BWVer = roipoly(ratio(:,:,analysisFrame),xTestVer, yTestVer);

    
    %%%%% Calculate the aveage within
     S1Mask(:,:,lauf) = img1(:,:,lauf) .*double(BWHor);
     Background3(lauf) = (sum(sum(S1Mask(:,:,lauf)))/BW_Box);

     S2Mask(:,:,lauf) = img2(:,:,lauf) .*double(BWHor);
     Background4(lauf) = (sum(sum(S2Mask(:,:,lauf)))/BW_Box);

     
    
    BW_Box_Out{subrun} = BW_Box;
    BWHor_Out{subrun} =BWHor;
   % BWVer_Out{subrun} = BWVer;
    x_Out{subrun} = x;
    y_Out{subrun} = y;
    
    end

end

lastFrameToConsider = pf;

Background1 =  mean(...
    reshape( ...
    img1(...
    (repmat(frapBkgIn, [1 1 lastFrameToConsider])) ... % 2d-mask in 3d: 2d mask replicated throughout time stack
    ), ...
    [], lastFrameToConsider));


Background2 =  mean(...
    reshape( ...
    img2(...
    (repmat(frapBkgIn, [1 1 lastFrameToConsider])) ... % 2d-mask in 3d: 2d mask replicated throughout time stack
    ), ...
    [], lastFrameToConsider));

for lauf= 1:pf

BackgroundA(lauf) = (Background1(lauf) + Background3(lauf)) ./2;
BackgroundB(lauf) = (Background2(lauf) + Background4(lauf)) ./2;

end

%(lauf)%Bkg1(:,:) = [Back(lauf)ground1' Backgr(lauf)ound3'];
%Bkg2((lauf):,:) = [Background2' Background4'];

%for lauf =1:pf
 %   imgCorr1(:,:,lauf) = img1(:,:,lauf) - uint16(BackgroundA(lauf));%%%% changed from uint16 to uint8
 %   imgCorr2(:,:,lauf) = img2(:,:,lauf) - uint16(BackgroundB(lauf));
%end