function [ r,theta,binedges] = polarPlotAnisoStat(aniso,orient,avg,varargin)

optargs.Parent=gca;
optargs.Nbins=20; %  Number of bins.
optargs.PlotType='XY';
optargs.ReferenceOrient=0; %Reference theta and anisotropy.
optargs=parsepropval(optargs,varargin{:});

% Bin edges and center.
binedges=(pi/180)*linspace(0,180,optargs.Nbins+1); % The last count by histc is the number of values that are exactly 180 degree.



% Counts
[~,indices]=histc(orient,binedges);

vectors=aniso.*exp(1i*2*orient);
meanVector=mean(vectors);
meanAniso=abs(meanVector);
meanOrient=mod(0.5*angle(meanVector),pi);

% Calculate bin-height by combining counts and weights.
binHeight=zeros(size(binedges));

for idBin=1:numel(binedges)
   % binHeight(idBin)=counts(idBin)*mean(weights(indices == idBin));
   binHeight(idBin)=sum(avg(indices == idBin)); % counts are accounted for by the sum of intensities.
end


% theta and r are variables used for plotting. 
theta=zeros(1,3*length(binedges));
r=zeros(1,3*length(binedges));

theta(1:3:end)=binedges;
theta(2:3:end)=binedges;
theta(3:3:end)=binedges;
r(3:3:end)=binHeight;
r(4:3:end)=binHeight(1:end-1);
switch(optargs.PlotType)
    case 'XY'
        plot(theta,r);
    case 'Polar'
        h=polar(theta,r,'m');
        hold on;
        XRange=get(gca,'XLim');
        XMax=XRange(2);
        hMean=polar([meanOrient meanOrient],[0 XMax],'r--');
        hRef=polar([optargs.ReferenceOrient optargs.ReferenceOrient],[0 XMax],'b--');
        set(hMean,'LineWidth',3);
        set(hRef,'LineWidth',3);
        hold off;
      
end
% 
% t=(pi/180)*(0:0.1*optargs.Interval:180)'; % t variable for polar plot. 
% r=zeros(size(t));
% 
% rho=zeros(size(binedges));
% for idt=1:numel(binedges)
%     rho(idt)=counts(idt)*mean(weights(indices==idt));
% end
% polar(optargs.Parent,binedges,rho,'k-');
% 


end

