%%%%%%% Structure Tensor implementation

getd = @(p)path(p,path);

getd('toolbox_signal/');
getd('toolbox_general/');

%n = 256;
name = 'test1.tif';


test = tiffread2('test1.tif');
test = test.data;

%[X,Y] = meshgrid(1:n,1:n);
%[XI,YI] = meshgrid(1:1/2:n+1/2,1:1/2:n+1/2); XI(:,end) = n; YI(end,:) = n;
%f1 = interp2(X,Y,f,XI,YI);

f1 = test;

G = grad(f1);
f = double(test);

n = 512;
%name = 'hibiscus';
%f = load_image(name,n);
f = rescale( sum(f,3) );
clf;
imageplot(f);

cconv = @(f,h)real(ifft2(fft2(f).*repmat(fft2(h),[1 1 size(f,3)])));

t = [0:n/2 -n/2+1:-1];
[X2,X1] = meshgrid(t,t);
normalize = @(h)h/sum(h(:));
h = @(sigma)normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );


blur = @(f,sigma)cconv(f,h(sigma));

options.order = 2;
nabla = @(f)grad(f,options);

tensorize = @(u)cat(3, u(:,:,1).^2, u(:,:,2).^2, u(:,:,1).*u(:,:,2));
rotate = @(T)cat(3, T(:,:,2), T(:,:,1), -T(:,:,3));

T = @(f,sigma)blur( tensorize( nabla(f) ), sigma);

options.sub = 8;
clf; sigma = 0.2;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);

clf; sigma = 4;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);
