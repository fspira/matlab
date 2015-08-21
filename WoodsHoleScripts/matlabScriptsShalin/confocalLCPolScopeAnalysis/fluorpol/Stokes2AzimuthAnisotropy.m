function [Azimuth, Anisotropy]=Stokes2AzimuthAnisotropy(S)
% Retrieve azimuth and anisotropy from stokes vector.

Azimuth=mod(0.5*atan2(S(3,:),S(2,:)),pi);
Anisotropy=sqrt(S(2,:).^2+S(3,:).^2)./S(1,:);
end
