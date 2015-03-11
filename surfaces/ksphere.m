function [fp,f] = ksphere(Am, Ap, Fm0, Fp0, Iy , Ix, mu, looporder)
%Updated 13 November 2014
%Computes the constant curvature surfaces in the 3-sphere using 
% the DPW (or generalized d'Alembert) method.  The computation is
%essentially the same as in the programme ksurf2, only a projection
% to S^3 is taken, instead of applying the Sym-formula.  The surface in
%S^3 is the output f, and the output fp is a stereographic projection
%to R^3, which is also plotted.

%The loop group construction, with the same setup used here, is
%described in 	arXiv:1301.5999 [math.DG]
%"Constant Gaussian curvature surfaces in the 3-sphere via loop groups"
%by David Brander, Jun-ichi Inoguchi, Shimpei Kobayashi. 

%All loops will be expanded up to order
%'looporder' or the maximum order of the other input data in the loop parameter,
%whichever is the greater.  Increasing 'looporder' will result in greater
%accuracy, but longer computation time.  Also plots the surface using surf.
%The loops described there
%have to be untwisted before entering them into this program.  This is
%achieved by: twisted  
%        [a(lambda),  b(lambda);  c(lambda),  d(lambda)]
%           ->
%        [a(sqrt(lambda)),  B_-(sqrt(lambda));  C_+(sqrt(lambda)), d(sqrt(lambda)]   
%  untwisted 
% where B_-(lambda):= lambda^{-1} b(lambda) and
%C_+(lambda):= lambda*c(lambda).


%Author: David Brander;  email: D.Brander@mat.dtu.dk.


% Input: All loops are untwisted. 
% Ap = @(t) [A_-n(t),... A_n(t)]
% where each A_i(t) is a 2 x 2 matrix valued function handle.
%This represents the potential (A_-n lambda^{-1} + A_1 lambda + ... + A_n
%lambda^n) dx. This corresponds to the "+" potential.
% Similarly, Am = @(t) [B_-m (t) .... B_0 (t)], 
% represents the minus potential (B_-m lambda^(-m) + ... + B_0] dy.

%"looporder" can be a positive integer "looporder=maxlooporder"
%or a pair "looporder=[maxlooporder, minlooporder]" 
%(Default is minlooporder = floor(maxlooporder/2)).  Loops
%are approximated as polynomials of order minlooporder at the initial 
%point up to maxlooporder at the points furthest from the initial point.
%Default for minlooporder is the 1 or the maximum order of Fp0 and Fm0.

% All loops are converted into large square matrices, such that every pair of rows is a copy of the previous
%pair, shifted to the right by 2 columns. Then matrix
% multiplication corresponds exactly to loop multiplication, and the LU
% decomposition corresponds to the birkhoff decomposition of loops.

%Fp_0 and Fm_0  represent the initial conditions, and are Laurent
%polynomial (untwisted) loops: Fp_0 = [F_-j , ...., F_j]
%   Fm_0 = [G_-k,  ....., G_k],   where F_i and G_i are constant 2x2
%   matrices.  Usually Fp_0 = Fm_0 = eye(2).
%Ix = [x_init, x_step, leftpoints, rightpoints], where x_init is the initial x-value, 
%x_step is the step size between the x-values at which the 
%solution should be evaluated, and the integration is done in both directions, leftpoints to the left and rightpoints to the right of x_0.
%If either of leftpoints or rightpoints is less than 3, then it will be
%replaced by 0.
%Iy is of a similar form.


%fp, (the projection to R^3) is of the form [X; Y; Z],
%where each of X,Y and Z are m x n matrices,
%m is the number of x gridpoints 1+Ix(3)+Ix(4), and n the number of y
%grid points 1+Iy(3)+Iy(4).
%X is the matrix of x-values at the grid
%points, Y the matrix of y values, and Z the matrix of z values. 
%This can be plotted with the command 
%   surf(f(1:m, :), f(m+1: 2*m, :),  f(2*m+1: 3*m, :)),
% (that is, surf(X,Y,Z) ) which is also done automatically.

%checkinput parses the input for consistency, converts the loops into
%square matrix form and sets the minimum and maximum looporders.
[xpoints,ypoints,bigAm, bigAp, bigFm0, bigFp0,minlooporder,maxlooporder]=checkinput(Am, Ap, Fm0, Fp0, Iy , Ix, looporder);
ordermatrix = setlooporders(Ix,Iy,minlooporder,maxlooporder); %the orders for the LU decomposition and loop evaluation
[Fplusinv,xerror]=invSintegrateA(bigAp, inv(bigFp0), Ix,2,1);
disp(['x-axis:  Left endpoint Error ', num2str(xerror(1),3),'. Right endpoint Error ', num2str(xerror(2),3),]);
FpDims=size(Fplusinv); 
bsize=FpDims(1);
[Fminus,yerror] = XSintegrateA(bigAm, bigFm0, Iy,2,1);
disp(['y-axis:  Left endpoint Error ', num2str(yerror(1),3),'. Right endpoint Error ', num2str(yerror(2),3),]);
Y=zeros(2*xpoints,2*ypoints);
ErrorMatrix=zeros(xpoints,ypoints);   
for j=1:xpoints
    factor1=looptranspose(Fplusinv(:,2*(j-1)+1:2*j),maxlooporder);
    for k=1:ypoints
        splitorder=ordermatrix(j,k);
           factor2=loop2SMat(Fminus(: ,(k-1)*bsize+2*maxlooporder+1:k*bsize-2*maxlooporder),splitorder,splitorder);
         Phi=factor1(:,2*(2*maxlooporder-splitorder)+1: 2*(2*maxlooporder+splitorder+1))*factor2;
         tosplit=loop2SMat(Phi,splitorder,splitorder);
       [~,right,~]=lu(tosplit,0.0000001);   
         F = Fminus(: ,(k-1)*bsize+2*(2*maxlooporder-splitorder)+1:k*bsize-2*(2*maxlooporder-splitorder))/right;        
       F1=rloopeval(F,1);
       Fmu=[1,mu;1/mu,1].*rloopeval(F,mu^2);
       f(2*j-1:2*j,2*k-1:2*k)=Fmu/F1;
       ErrorMatrix(j,k)=(norm(f(2*j-1:2*j,2*k-1:2*k)'*f(2*j-1:2*j,2*k-1:2*k)-eye(2)))/sqrt(2);
    end
    if mod(j,20)==1
            disp(['Row ', num2str(j),  ' Max Error ', num2str(max(ErrorMatrix(j:j,:)),2),'. Errors: ', num2str(ErrorMatrix(j:j,:),2)]);
    end      
end
disp(' ');
disp(['Max error:', num2str(max(max(ErrorMatrix)),2), '. Mean error: ', num2str(mean(mean(ErrorMatrix)),2), '. Corner Errors: ']);
disp(num2str([ErrorMatrix(1,1), ErrorMatrix(1,ypoints);ErrorMatrix(xpoints,1), ErrorMatrix(xpoints,ypoints)],2));
fp = projplotdata(f);
end
    

     %%%%%%%%%%%%%%%%%%%%%%%%%%%% sub functions %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = rloopeval(X,r)
% Evaluates a loop at lambda = r.
% Input: X = [x_-n .... x_0 x_1 .... x_n], where each x_i is 2x2
%  r is a number.
%Output: a sparse 2x2 matrix, the sum of x_i r^i
n = (length(X)-2)/4; % the order of X in lambda
Y= full(X(1:2, 2*n+1:2*n+2));
for count = 1 : n 
    Y = Y +  r^count*X(1:2, 2*n+1+2*count : 2*n+2+ 2*count) +  r^(-count)*X(1:2, 2*n+1-2*count : 2*n + 2 - 2*count);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = projplotdata(f)
%
fheight=length(f(:,1:1));
xpoints=fheight/2;

fwidth=length(f(1:1,:));
ypoints=fwidth/2;

Y=zeros(3*xpoints, ypoints);
for count=1:xpoints
    for j=1:ypoints
        Y(count,j) = -imag(f(2*count-1, 2*j-1))/(1+real(f(2*count-1, 2*j-1)));    
        Y(xpoints + count, j) = -real(f(2*count-1, 2*j))/(1+real(f(2*count-1, 2*j-1))); 
        Y(2*xpoints+count, j) =  -imag(f(2*count-1, 2*j))/(1+real(f(2*count-1, 2*j-1)));         
    end
end
        
surf(Y(1:xpoints, :), Y(xpoints+1: 2*xpoints, :),  Y(2*xpoints+1: 3*xpoints, :),'EdgeAlpha',0.3,'FaceAlpha',1,'FaceLighting','gouraud','EdgeLighting','gouraud','EdgeColor','white','FaceColor',[0.05 0.4 0.8])
daspect([1 1 1])
grid off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [Finv,Error] = invSintegrateA(A, F0, Ix,m,lambda)
%The integration along the y axis.
%solves the equation dG/dt = -A(t) G

%The values of Finv are returned as loops in column form.
% Estimates the "error" at the endpoints  by checking if the determinant is
% 1.

tempIy3=Ix(3); tempIy4=Ix(4);
if Ix(3)<2, tempIy3=2; end
if Ix(4)<2, tempIy4=2; end
stepy= Ix(2);   intervalyr=zeros(1,tempIy4+1); intervalyl=zeros(1,tempIy3+1);
for j = 1 : tempIy4+1
    intervalyr(j)=  Ix(1) + (j-1)*stepy;
end
for j = 1 : tempIy3+1
    intervalyl(j)=  Ix(1) - (j-1)*stepy;
end
bsize=length(F0);
n =(bsize/m-1)/4; % the order of A in lambda
initF=reshape(F0(:, m*2*n+1:m*2*n+m),m*bsize,1);
testfn = @(t,y) reshape((-A(t))*sparse(reshape(y,bsize,m)), m*bsize,1);
if Ix(3)>0
  [~,soll] = ode45(testfn ,intervalyl,initF);
  G1=sparse(reshape(flipud(soll(2:Ix(3)+1, :)).', bsize,m*(Ix(3))));
else
 G1=[];
end
if Ix(4)>0
  [~,solr] = ode45(testfn ,intervalyr,initF);
  G2=sparse(reshape(solr(2:Ix(4)+1, :).', bsize,m*(Ix(4))));
else
    G2=[];
end
Finv=[G1, F0(:, m*2*n+1:m*2*n+m), G2];
Error=zeros(2,1);
Error(1)=abs(1-det(rloopeval(looptranspose(Finv(:,1:m),n),lambda)));
Error(2)=abs(1-det(rloopeval(looptranspose(Finv(:,m*(Ix(3)+Ix(4))+1:m*(Ix(3)+Ix(4)+1)),n),lambda)));
end

   
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [F,Error] = XSintegrateA(A, F0, Ix,m,lambda)
%solves the equation dF/dt = F A, with initial condition F(t0) = F0, with t
%in the specified interval.  
% Estimates the "error" at the endpoints  by checking if the determinant is
% 1.

tempIx3=Ix(3); tempIx4=Ix(4);
if Ix(3)<2, tempIx3=2; end
if Ix(4)<2, tempIx4=2; end

stepy= Ix(2);  intervalyl=zeros(1, tempIx3+1); intervalyr=zeros(1,tempIx4+1); 
for count = 1 : tempIx4+1
    intervalyr(count)=  Ix(1) + (count-1)*stepy;
end
for count=1: tempIx3+1
    intervalyl(count)=  Ix(1) - (count-1)*stepy;
end
bsize=length(F0);
n =(bsize/m-1)/4; % the order of X in lambda
initF=reshape(F0(m*2*n+1:m*2*n+m,:),m*bsize,1);
y=zeros(m*bsize,1);
testfn = @(t,y) reshape(reshape(y,m,bsize)*A(t), m*bsize,1);
[~,sol1] = ode45(testfn ,intervalyl,initF); 
[~,sol2] = ode45(testfn ,intervalyr,initF);  
F=[reshape(flipud(sol1(2:Ix(3)+1, :)).', m,bsize*(Ix(3))),   reshape(sol2(1:Ix(4)+1, :).', m,bsize*(Ix(4)+1))];
Error=zeros(2,1);
Error(1)=abs(1-det(rloopeval(F(:,1:bsize),lambda)));
Error(2)=abs(1-det(rloopeval(F(:,bsize*(Ix(3)+Ix(4))+1:bsize*(Ix(3)+Ix(4)+1)),lambda)));
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V=looptranspose(U,d)
   %input: U is a loop given as a block colum matrix
    % U= [U_(n);  ...  ; U_(-n)]
    % output:  V is in row form V= [0,...,0,V_(-n+d),  ... V_(n-d),0,...0]
    Udims=size(U);
    len=Udims(1);
    m=Udims(2);
    K=len/m;  %the number of coefficients of U
    V=sparse(zeros(m,m*K));
      for j=d:K-1-d
         V(:,m*j+1:m*j+m)=U(len-m*(j+1)+1:len-m*(j+1)+m,:);
      end       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xmat = loop2matrix(X,n,k)
%Converts a loop into a big square matrix 
%Input:  a matrix X=[x_-nold .... x_0 x_1 .... x_nold] representing a Laurent
%polynomial loop of order nold
%The input loop does not need to be the same order as the output loop

%x_-nold lambda^-nold + ... + x_0 + x_nold lambda^nold, where each x_i is an mxm matrix 
%Output:  Xmat is a sparse m(4n+1) x m(4n+1) square matrix, which looks something
%like  [x_0 ... x_n 0 0 ... 0;
%       x_-1 x_0 ....x_n 0 ...;
%       x_-2 x_-1 x_0 .......; ..... etc
% Multiplying  Xmat with Ymat corresponds to multiplying the loops X and Y.
Xdims=size(X);
m=Xdims(1);
nold =(Xdims(2)/m-1)/2; % the order of the input loop X in lambda
k=min(nold,k);
Row=X(1:m, m*(nold-k)+1: m*(nold+k+1));
Xmat=sparse(zeros(m*(4*n+1)));
for j=0:n
    Xmat(m*n+m*j+1:m*n+m*j+m, m*(n-k+j)+1:m*(n-k+j+2*k+1))=Row;
end
Xmat(1:m*n, 1:m*(2*n+1))=Xmat(m*n+1:2*m*n,m*n+1:m*3*n+m);
Xmat(m*(2*n+1)+1:m*(3*n+1), m*(n)+1:m*(4*n+1))=Xmat(m*n+m+1:m*(2*n+1),1:m*(3*n+1));
Xmat(m*(3*n+1)+1:m*(4*n+1), m*(2*n+1)+1:(4*n+1)*m)=Xmat(m*n+1:2*m*n,1:2*m*n);
end           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xmat = loop2SMat(X,n,k)
%Converts a loop into a big square matrix 
%Input:  a matrix X=[x_-nold .... x_0 x_1 .... x_nold] representing a Laurent
%polynomial loop of order nold
%The input loop does not need to be the same order as the output loop

%x_-nold lambda^-nold + ... + x_0 + x_nold lambda^nold, where each x_i is an mxm matrix 
%Output:  Xmat is a sparse m(4n+1) x m(4n+1) square matrix, which looks something
%like  [x_0 ... x_n 0 0 ... 0;
%       x_-1 x_0 ....x_n 0 ...;
%       x_-2 x_-1 x_0 .......; ..... etc
% Multiplying  Xmat with Ymat corresponds to multiplying the loops X and Y.
Xdims=size(X);
m=Xdims(1);
nold =(Xdims(2)/m-1)/2; % the order of the input loop X in lambda
k=min(nold,k);
Xmat=zeros(m*(2*n+1));
Xmat(m*n+1:m*(n+1),m*(n-k)+1:m*(n+1+k))=full(X(1:m, m*(nold-k)+1: m*(nold+k+1)));
for j=0:n-1
    Xmat(m*j+1:m*j+m,1:m*(n+1+j))=Xmat(m*n+1:m*(n+1),m*(n-j)+1:m*(2*n+1));
    Xmat(m*(2*n-j)+1:m*(2*n-j+1),m*(n-j)+1:m*(2*n+1))=Xmat(m*n+1:m*(n+1),1:m*(n+1+j));
end
Xmat=sparse(Xmat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%

    function orders = setlooporders(Ix,Iy,minlooporder,maxlooporder)
 %creates a matrix of looporders to be used in the LU decomposition
 %to speed things up.
 
orders=zeros(Ix(3)+Ix(4)+1,Iy(3)+Iy(4)+1);
center=[Ix(3)+1,Iy(3)+1];            %the coordinates of the center of the matrix
cornerdistances=[norm(center-[1,1]), norm(center-[1,Iy(3)+Iy(4)+1]), norm(center-[Ix(3)+Ix(4)+1,1]),  norm(center-[Ix(3)+Ix(4)+1,Iy(3)+Iy(4)+1])];
for m=1:(Ix(3)+Ix(4)+1)
   for k=1:Iy(3)+Iy(4)+1
        orders(m,k)=min(ceil(minlooporder+(maxlooporder-minlooporder)*norm(center-[m,k])/max(cornerdistances)),maxlooporder);
   end
end
     end 
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xpoints,ypoints,bigAm, bigAp, bigFm0, bigFp0,minlooporder,maxlooporder]=checkinput(Am, Ap, Fm0, Fp0, Iy , Ix, looporder)
%Checks the input data for consistency with expected input.
%Converts the intial data, in terms of simple untwisted loops, into 
%untwisted loops represented by big square matrices, expanded up to order
%'looporder' or the maximum order of the other data in the loop parameter,
%whichever is the greater.

% Output:  Each of bigAp and bigAm is a handle for a function of the form  = @(t) M(t)
% where M is a  2(4n+1) x 2(4n+1) square matrix valued function of t
% which looks something like 
%      [x_0 ... x_n 0 0 ... 0;
%       x_-1 x_0 ....x_n 0 ...;
%       x_-2 x_-1 x_0 .......; ..... etc,  each x_i a 2x2 matrix

% bigFp0 and bigFm0 are square matrices of initial values, of the same size 
%and similar structure to A.
assert((Ix(3)>=0)&(Ix(4)>=0)&(Ix(3)+Ix(4)>0),'Ix(3) and Ix(4) must be non-negative, and at least one positive');
assert((Iy(3)>=0)&(Iy(4)>=0)&(Iy(3)+Iy(4)>0),'Iy(3) and Iy(4) must be non-negative, and at least one positive');
xpoints=1+Ix(3)+Ix(4); ypoints=1+Iy(3)+Iy(4);

Fm0Dims=size(Fm0); 
Fp0Dims=size(Fp0);
AmDims=size(Am(Iy(1)));
ApDims=size(Ap(Ix(1)));
m=Fm0Dims(1);
assert((Fm0Dims(1)==2)&(Fp0Dims(1)==2)&(AmDims(1)==2)&(ApDims(1)==2),'All input loops must be matrices with two rows');

Fm0order = (Fm0Dims(2)-m)/(2*m);  Fp0order = (Fp0Dims(2)-m)/(2*m);
AmOrder = (AmDims(2)-m)/(2*m);  ApOrder = (ApDims(2)-m)/(2*m);
assert(norm(floor([Fm0order, Fp0order, AmOrder,ApOrder])-[Fm0order, Fp0order, AmOrder,ApOrder])==0,...
    'Dimensions of input matrices and function handles are not consistent for loops');
loDims=size(looporder);
assert((loDims(2)<=2)&(norm(floor(looporder)-looporder)==0)&(norm(looporder-abs(looporder))==0),...
    'looporder must be a non-negative integer or a 1x2 non-negative integer matrix ');
if loDims(2)==2
    maxlooporder=max([looporder(1),Fm0order,Fp0order,AmOrder,ApOrder]);
    minlooporder=looporder(2);
else
  maxlooporder=max([looporder,Fm0order,Fp0order,AmOrder,ApOrder]); 
  minlooporder=max([floor(maxlooporder/2),Fm0order,Fp0order]);  
end
%convert the loops into square matrix form:
bigAm=@(t)loop2matrix(Am(t),maxlooporder,AmOrder);
bigAp=@(t)loop2matrix(Ap(t),maxlooporder,ApOrder);
bigFm0= loop2matrix(Fm0,maxlooporder,Fm0order);
bigFp0= loop2matrix(Fp0,maxlooporder,Fp0order);

 end


