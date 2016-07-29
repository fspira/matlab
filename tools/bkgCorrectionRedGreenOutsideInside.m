function [BackgroundA,BackgroundB,Background1, Background2, Background3, Background4,imgCorr1,imgCorr2] = bkgCorrectionRedGreenOutsideInside(img1,img2)

%img1 = greenStack;

imgNorm = normalizedImage(img2(:,:,1));

[mf nf pf] = size(img2);

fh = figure(1);
 title('Mark the Background')
frapBkgIn = roipoly(imgNorm(:,:,1));
close(fh);



fh = figure(1);
 title('Mark the Background')
frapBkgOut = roipoly(imgNorm(:,:,1));
close(fh);

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


Background3 =  mean(...
    reshape( ...
    img1(...
    (repmat(frapBkgOut, [1 1 lastFrameToConsider])) ... % 2d-mask in 3d: 2d mask replicated throughout time stack
    ), ...
    [], lastFrameToConsider));


Background4 =  mean(...
    reshape( ...
    img2(...
    (repmat(frapBkgOut, [1 1 lastFrameToConsider])) ... % 2d-mask in 3d: 2d mask replicated throughout time stack
    ), ...
    [], lastFrameToConsider));


BackgroundA = (Background1 + Background3) ./2;
BackgroundB = (Background2 + Background4) ./2;

%Bkg1(:,:) = [Background1' Background3'];
Bkg2(:,:) = [Background2' Background4'];

for lauf =1:pf
    imgCorr1(:,:,lauf) = img1(:,:,lauf) - uint16(BackgroundA(lauf));%%%% changed from uint16 to uint8
    imgCorr2(:,:,lauf) = img2(:,:,lauf) - uint16(BackgroundB(lauf));
end