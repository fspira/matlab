function normFactors = computeConfocalNormalizations( isopath,Wavelength, NA,PixSize )

if(isempty(isopath))
 isopath=uigetdir('Select folder tha contains isotropic calibration');
 if(isempty(isopath))
     normFactors=[1 1 1 1];
     return;
 end
end

[I0Iso,I45Iso,I90Iso,I135Iso,~,blackiso]=readconfocalPolData([isopath '/P001/']);
sigma=0.2*Wavelength/NA;
sigmaPix=sigma/PixSize;
FiltGauss=fspecial('gaussian',round(7*sigmaPix),sigmaPix);

I0IsoFilt=imfilter(I0Iso,FiltGauss,'replicate','same');
I45IsoFilt=imfilter(I45Iso,FiltGauss,'replicate','same');
I90IsoFilt=imfilter(I90Iso,FiltGauss,'replicate','same');
I135IsoFilt=imfilter(I135Iso,FiltGauss,'replicate','same');

eq45=I0IsoFilt./I45IsoFilt;
eq90=I0IsoFilt./I90IsoFilt;
eq135=I0IsoFilt./I135IsoFilt;

normFactors={ones(size(eq45)), eq45, eq90, eq135};

[OrientationBeforeNormalization, AnisotropyBeforeNormalization, IntensityBeforeNormalization]=...
    ComputeFluorAnisotropy(I0Iso,I45Iso,I90Iso,I135Iso,...
    'anisotropy','BlackLevel',0,'normFactors',[1 1 1 1]);

[OrientationAfterNormalization, AnisotropyAfterNormalization, IntensityAfterNormalization]=...
    ComputeFluorAnisotropy(I0Iso,I45Iso,I90Iso,I135Iso,...
    'anisotropy','BlackLevel',0,'normFactors',normFactors);

hwideIso=togglefig('isotropic slide',1); colormap gray;
set(hwideIso,'Position',[100 100 2500 1000],'defaultaxesfontsize',15);

ha=imagecat(OrientationBeforeNormalization,OrientationAfterNormalization, OrientationBeforeNormalization,OrientationAfterNormalization, AnisotropyBeforeNormalization, AnisotropyAfterNormalization, IntensityBeforeNormalization,IntensityAfterNormalization,...
    'equal','colorbar','off');

[countsIso,levelsIso]=hist((180/pi)*OrientationBeforeNormalization(:),0:0.5:180);
axes(ha(3)); stem(levelsIso,countsIso);  title('Histogram of orientaiton before normalization');
axis tight; xlim([0 180]); xlabel('Orientation');
 
[countsIsoEq,levelsIsoEq]=hist((180/pi)*OrientationAfterNormalization(:),0:0.5:180);
axes(ha(4)); stem(levelsIsoEq,countsIsoEq);  title('Histogram of orientaiton after normalization');
axis tight; xlim([0 180]); xlabel('Orientation');


end

