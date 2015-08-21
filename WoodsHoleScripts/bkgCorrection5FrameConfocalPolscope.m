function [imgCorr1,imgCorr2,imgCorr3,imgCorr4,imgCorr5] =  bkgCorrection5FrameConfocalPolscope(img1,img2,img3,img4,img5)
%img1 = I0File(:,:,1);
imgNorm = normalizedImage(img1(:,:,1));

[mf nf pf] = size(img1);

fh = figure(1);
 title('Mark the Background')
frapBkg = roipoly(imgNorm(:,:,1));
close(fh);

pixelInArea = sum(sum(frapBkg));

if class(img1) == 'uint32'

bkgTmp  = uint32(frapBkg) .* img1(:,:,1);
Background1= sum(sum(bkgTmp))/ pixelInArea


lastFrameToConsider = pf;

for lauf =1:pf
    imgCorr1(:,:,lauf) = uint32(img1(:,:,lauf)) - round(Background1);
    imgCorr2(:,:,lauf) = uint32(img2(:,:,lauf)) - round(Background1);
    imgCorr3(:,:,lauf) = uint32(img3(:,:,lauf)) - round(Background1);
    imgCorr4(:,:,lauf) = uint32(img4(:,:,lauf)) - round(Background1);
    imgCorr5(:,:,lauf) = uint32(img5(:,:,lauf)) - round(Background1);
end


elseif class(img1) == 'uint16'
    
    
bkgTmp  = uint16(frapBkg) .* img1(:,:,1);
Background1= sum(sum(bkgTmp))/ pixelInArea


lastFrameToConsider = pf;

for lauf =1:pf
    imgCorr1(:,:,lauf) = uint16(img1(:,:,lauf)) - round(Background1);
    imgCorr2(:,:,lauf) = uint16(img2(:,:,lauf)) - round(Background1);
    imgCorr3(:,:,lauf) = uint16(img3(:,:,lauf)) - round(Background1);
    imgCorr4(:,:,lauf) = uint16(img4(:,:,lauf)) - round(Background1);
    imgCorr5(:,:,lauf) = uint16(img5(:,:,lauf)) - round(Background1);
end

elseif class(img1) == 'uint8'
    
    bkgTmp  = uint8(frapBkg) .* img1(:,:,1);
Background1= sum(sum(bkgTmp))/ pixelInArea


lastFrameToConsider = pf;

for lauf =1:pf
    imgCorr1(:,:,lauf) = uint8(img1(:,:,lauf)) - round(Background1);
    imgCorr2(:,:,lauf) = uint8(img2(:,:,lauf)) - round(Background1);
    imgCorr3(:,:,lauf) = uint8(img3(:,:,lauf)) - round(Background1);
    imgCorr4(:,:,lauf) = uint8(img4(:,:,lauf)) - round(Background1);
    imgCorr5(:,:,lauf) = uint8(img5(:,:,lauf)) - round(Background1);
end

end