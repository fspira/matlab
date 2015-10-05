function [orient, aniso, avg]=ComputeFluorAnisotropy_1(I0,I45,I90,I135,process,varargin)
% [orient,aniso, avg]=ComputeFluorAnisotropy(I0,I45,I90,I13,process,instrument matrix) computes magnitude
% and azimuth of fluorescence anisotropy.

% Author and Copyright: Shalin Mehta, HFSP Postdoctoral Fellow
%                           Marine Biological Laboratory, Woods Hole, MA
%                           http://www.mshalin.com
% 
% License: Restricted academic use. 
% This software is free to use, but only with explicit permission by the
% author and only for academic purpose. This software constitutes a part of
% unpublished method for measuring orientation of single and ensemble of
% fluorophores.
I0=double(I0);
I45=double(I45);
I90=double(I90);
I135=double(I135);
sizeOrig=size(I0);

optargs.ItoSMatrix=[0.5 0.5 0.5 0.5; 1 0 -1 0; 0 1 0 -1];
optargs.BlackLevel=0;
optargs.normFactors=[1 1 1 1]; % Factors for normalizing detection and/or excitation throughput.
optargs.anisoCeiling=1;
optargs.anisoFloor=0;
optargs.OrientReference=0; % Orientation reference angle in degrees.
optargs=parsepropval(optargs,varargin{:});

% Create variables for the sake of brevity.
ItoSMatrix=optargs.ItoSMatrix;
BlackLevel=optargs.BlackLevel;
normFactors=optargs.normFactors;
%%%%

% Normalize
if(isnumeric(normFactors))
    I0column=normFactors(1)*(I0(:)-BlackLevel);
    I45column=normFactors(2)*(I45(:)-BlackLevel);
    I90column=normFactors(3)*(I90(:)-BlackLevel);
    I135column=normFactors(4)*(I135(:) - BlackLevel);
elseif(iscell(normFactors))
    if(~ (ismatrix(normFactors{1}) && ismatrix(normFactors{2})  && ismatrix(normFactors{3})  && ismatrix(normFactors{4})) )
        error('Normalization factors expected to be 2D images of the same size as data.');
    end
    dim1=size(I0,1); dim2=size(I0,2); dim3=size(I0,3); dim4=size(I0,4); dim5=size(I0,5); % Dimensions of hyperstack.
    
    normFactorI0=imresize(normFactors{1},[dim1 dim2],'bilinear');
    normFactorI0=repmat(normFactorI0,1,1,dim3,dim4,dim5);
    
    normFactorI45=imresize(normFactors{2},[dim1 dim2],'bilinear');
    normFactorI45=repmat(normFactorI45,1,1,dim3,dim4,dim5);
    
    normFactorI90=imresize(normFactors{3},[dim1 dim2],'bilinear');
    normFactorI90=repmat(normFactorI90,1,1,dim3,dim4,dim5);
    
    normFactorI135=imresize(normFactors{4},[dim1 dim2],'bilinear');
    normFactorI135=repmat(normFactorI135,1,1,dim3,dim4,dim5);
    
    I0column=normFactorI0(:).*(I0(:)-BlackLevel);
    I45column=normFactorI45(:).*(I45(:)-BlackLevel);
    I90column=normFactorI90(:).*(I90(:)-BlackLevel);
    I135column=normFactorI135(:).*(I135(:) - BlackLevel);    
else
    error('Normalization factors should either be scalar factors or the cell aray of normalizing images.');
end

S=ItoSMatrix*[I0column I45column I90column I135column]';

orient=mod(0.5*atan2(S(3,:),S(2,:)),pi);
avg=0.5*S(1,:);
switch(process)
    case('anisotropy')
        % anisotropy=(Imax-Imin)/(Imax+Imin). Value ranges from 0 to 1.
     aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
     
    case('ratio')
        % ratio=Imax/Imin. Ranges from 1 to inf.
        aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
         aniso=(1+aniso)./(1-aniso); %Then Ratio
         aniso(aniso<1)=NaN;
         
    case('difference')
        % difference=Imax-Imin.
        aniso=sqrt(S(3,:).^2 + S(2,:).^2);
end
        
orient=reshape(orient,sizeOrig)-optargs.OrientReference;
avg=reshape(avg,sizeOrig);
aniso=reshape(aniso,sizeOrig);

% Clip the anisotropy above ceiling.
aniso(aniso>optargs.anisoCeiling)=optargs.anisoCeiling;

end