function  [D, dNegInv, dSum] = doVisualizeSPComponents(D,tifFilename)

[m n p] = size(D)
for lauf = 1:p
    minimumValue = min(min(D(:,:,lauf)));
    dTmp = D(:,:,lauf) > minimumValue & D(:,:,lauf) <0;
    dNeg = D(:,:,lauf) .* dTmp;
    dNegInv(:,:,lauf) = -dNeg;
end

%%%%%% Remove negative Values from D
tmp = D >0;
D = D .* tmp;


%%%%%% replace inf with zero
D(~isfinite(D)) = 0;
dNegInv(~isfinite(dNegInv)) = 0;


dNorm = normalizedImage3D(D);
dNegNorm = normalizedImage3D(dNegInv);

for lauf = 1:p
    dNegRange(:,:,lauf) = (dNegInv(:,:,lauf)./max(max(dNegInv(:,:,lauf)))).*127;
    dRange(:,:,lauf) = ((D(:,:,lauf)./ max(max(D(:,:,lauf)))).*126)+128;
    tmp = dRange(:,:,lauf) >129;
    dRange(:,:,lauf) = dRange(:,:,lauf) .*tmp;
end

%ratio = dNorm./dNegNorm;
dSum = dNegRange +dRange;
imshow(dSum(:,:,1),[]);
dSum = uint8(round(dSum));
colormap(jet)
figure(1)
imshow(dNorm(:,:,1),[])
figure(2)
imshow(dNegNorm(:,:,1),[])

%imwrite(dSum,[tifFilename,'_Aniso_FalseColor.tif'])
tiffwrite_mat(dSum, [tifFilename,'_Aniso_FalseColor.tif']);

%colormap(jet)
close all

tiffwrite_mat(dNegInv, [tifFilename,'_P.tif']);
tiffwrite_mat(D, [tifFilename,'_S.tif']);