function A_out=TV_A3(A_p,Y,X,rho,opts5)
if isfield(opts5,'mu');         mu = opts5.mu;                   else  mu = 10;                   end
if isfield(opts5,'beta');       beta= opts5.beta;             else  beta = 1000;               end
if isfield(opts5,'F_it');       maxit = opts5.F_it;                else  maxit = 15;                  end
m=size(A_p,1);
n=size(A_p,2);

[U1,S1,~]=svd(X*X');
Sig1=mu*diag(S1);

diaga=ones(m,1);diagb=ones(m-1,1);
D=diag(-diaga)+diag(diagb,1);
D(end,1)=1;

d=D(:,1);
deig=fft(d);
Sig2=beta*(abs(deig).^2);

Sig=repmat(Sig1',m,1)+repmat(Sig2,1,n)+rho;
Sig=1./Sig;

Theta=zeros(m,n);
W=zeros(m,n);

for i=1:maxit
    %% A subproblem    
    M=mu*Y*X'+beta*D'*(W-Theta)+rho*A_p;
    temp=Sig.*(fft(M)*U1);
    A_out=real(ifft(temp))*U1';
    %% W subproblem
    Thresh=1/beta;
    W=wthresh(D*A_out+Theta,'s',Thresh);
    %% updating Theta
    Theta=Theta+D*A_out-W;
end
end