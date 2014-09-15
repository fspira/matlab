n = 512;
name = 'fingerprint';
f = rescale(load_image(name,n));
clf;
imageplot(f);


test = tiffread2('test1.tif');
test = test.data;

%[X,Y] = meshgrid(1:n,1:n);
%[XI,YI] = meshgrid(1:1/2:n+1/2,1:1/2:n+1/2); XI(:,end) = n; YI(end,:) = n;
%f1 = interp2(X,Y,f,XI,YI);

f1 = test;


[X,Y] = meshgrid(1:n,1:n);
[XI,YI] = meshgrid(1:1/2:n+1/2,1:1/2:n+1/2); XI(:,end) = n; YI(end,:) = n;
f1 = interp2(X,Y,f,XI,YI);

G = grad(f1);

T0 = zeros(n*2,n*2,2,2);
T0(:,:,1,1) = G(:,:,2).^2;
T0(:,:,2,2) = G(:,:,1).^2;
T0(:,:,1,2) = -G(:,:,1).*G(:,:,2);
T0(:,:,2,1) = -G(:,:,1).*G(:,:,2);

sigma = 5;
T = perform_blurring(T0, sigma);

T = T(1:2:end,1:2:end,:,:);

options.sub = 8;
clf;
plot_tensor_field(T, [], options);

[e1,e2,lambda1,lambda2] = perform_tensor_decomp(T);

E = sqrt(lambda1+lambda2);
A = sqrt(lambda1./lambda2);

clf;
imageplot({E A}, {'Energy', 'Anisotropy'});

rho = .1;

T = perform_tensor_decomp(e1,e2,ones(n),ones(n)*rho);

TensorMult = @(T,u)cat(3,  T(:,:,1,1).*u(:,:,1) + T(:,:,1,2).*u(:,:,2), ...
                            T(:,:,2,1).*u(:,:,1) + T(:,:,2,2).*u(:,:,2) );
                        
                        dt = .2;
                        
                        f1 = f;
                        
                        f1 = f1 + dt * div( TensorMult(T, grad(f1) ) );
                    
                        I = find( E<=mean(E(:))*1.1 );
T1 = T;
T1([I; I+n^2; I+2*n^2; I+3*n^2]) = 0;
options.sub = 4;
clf;
plot_tensor_field(T1, [], options);

t = 20;
dt = .2;
niter=  round(t/dt);
kdisp = round(linspace(0,niter,5)); kdisp(1) = [];
k =1;
f1 = f; % rand(n);
for i=1:niter
    f1 = f1 + dt * div( TensorMult(T, grad(f1,options) ) );
    if i==kdisp(k)
        subplot(2,2,k);
        imageplot(f1);
        k = k+1;
    end
end

t = 20;
dt = .2;
niter=  round(t/dt);
kdisp = round(linspace(0,niter,5)); kdisp(1) = [];
f1 = rand(n);
v = []; k =1;
for i=1:niter
    f1 = f1 + dt * div( TensorMult(T, grad(f1,options) ) );
    f1 = perform_hist_eq(f1,'linear');
    if i==kdisp(k)
        subplot(2,2,k);
        imageplot(f1);
        k = k+1;
    end
end

[X,Y] = meshgrid(1:n,1:n);
[XI,YI] = meshgrid(1:1/2:n+1/2,1:1/2:n+1/2); XI(:,end) = n; YI(end,:) = n;
f1 = interp2(X,Y,f,XI,YI);
% Compute the gradient field.
G = grad(f1);
% Compute the rank-1 tensor field associated to the gradient rotated by
% 90°.
T0 = zeros(n*2,n*2,2,2);
T0(:,:,1,1) = G(:,:,2).^2;
T0(:,:,2,2) = G(:,:,1).^2;
T0(:,:,1,2) = -G(:,:,1).*G(:,:,2);
T0(:,:,2,1) = -G(:,:,1).*G(:,:,2);
% Blur the tensor field.
sigma = 20;
T = perform_blurring(T0, sigma);
% Sub-sample it.
T = T(1:2:end,1:2:end,:,:);
options.sub = 4;
clf;
plot_tensor_field(T, [], options);