function [imgCorr1,imgCorr2] = bkgCorrectionRedGreen(img1,img2)

imgNorm = normalizedImage(img1(:,:,1));

[mf nf pf] = size(img1);

fh = figure(1);
 title('Mark the Background')
frapBkg = roipoly(imgNorm(:,:,1));
close(fh);

lastFrameToConsider = pf;

Background1 =  mean(...
    reshape( ...
    img1(...
    (repmat(frapBkg, [1 1 lastFrameToConsider])) ... % 2d-mask in 3d: 2d mask replicated throughout time stack
    ), ...
    [], lastFrameToConsider));


Background2 =  mean(...
    reshape( ...
    img2(...
    (repmat(frapBkg, [1 1 lastFrameToConsider])) ... % 2d-mask in 3d: 2d mask replicated throughout time stack
    ), ...
    [], lastFrameToConsider));

for lauf =1:pf
    imgCorr1(:,:,lauf) = img1(:,:,lauf) - uint16(Background1(lauf));
    imgCorr2(:,:,lauf) = img2(:,:,lauf) - uint16(Background2(lauf));
end