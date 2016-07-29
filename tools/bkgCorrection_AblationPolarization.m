function [imgCorr1,imgCorr2 frapBkg]  = bkgCorrection_AblationPolarization(img, img1)

imgNorm = normalizedImage(img(:,:,2));
imgNorm1 = normalizedImage(img1(:,:,2));

[mf nf pf] = size(img);

fh = figure(1);
 title('Mark the Background')
frapBkg = roipoly(imgNorm(:,:,1));
close(fh);

lastFrameToConsider = pf;

%for funcRun = 1:lastFrameToConsider
    
    Background1=  sum(sum((img(:,:,2) .* uint8(frapBkg))))/sum(sum(frapBkg))

    Background2 = sum(sum((img1(:,:,2) .* uint8(frapBkg))))/sum(sum(frapBkg))

%end



for lauf =1:pf
    imgCorr1(:,:,lauf) = img(:,:,lauf) - double(Background1);
    imgCorr2(:,:,lauf) = img1(:,:,lauf) - double(Background2);
end