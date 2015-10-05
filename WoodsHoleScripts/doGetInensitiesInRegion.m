function [I0Values,I45Values,I90Values,I135Values,IBValues] = doGetInensitiesInRegion(I0File,I45File,I90File,I135File,IBFile)


imgNorm = normalizedImage(I0File(:,:,1));


fh = figure(1);
 title('Mark the Background')
frapBkg = roipoly(imgNorm(:,:,1));
close(fh);

pixelInArea = sum(sum(frapBkg));


bkgTmp  = uint8(frapBkg) .* I0File(:,:,1);
I0Values = sum(sum(bkgTmp))/ pixelInArea



bkgTmp  = uint8(frapBkg).* I45File(:,:,1);
I45Values= sum(sum(bkgTmp))/ pixelInArea



bkgTmp  = uint8(frapBkg) .* I90File(:,:,1);
I90Values= sum(sum(bkgTmp))/ pixelInArea


bkgTmp  = uint8(frapBkg) .* I135File(:,:,1);
I135Values= sum(sum(bkgTmp))/ pixelInArea


bkgTmp  = uint8(frapBkg) .* IBFile(:,:,1);
IBValues= sum(sum(bkgTmp))/ pixelInArea




end
