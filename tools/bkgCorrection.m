function imgCorr = bkgCorrection(img)

imgNorm = normalizedImage(img(:,:,1));

[mf nf pf] = size(img);

fh = figure(1);
 title('Mark the Background')
frapBkg = roipoly(imgNorm(:,:,1));
close(fh);

lastFrameToConsider = pf;

Background =  mean(...
    reshape( ...
    img(...
    (repmat(frapBkg, [1 1 lastFrameToConsider])) ... % 2d-mask in 3d: 2d mask replicated throughout time stack
    ), ...
    [], lastFrameToConsider));

for lauf =1:pf
    imgCorr(:,:,lauf) = img(:,:,lauf) - double(Background(lauf));
end