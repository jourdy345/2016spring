%Simulate and save the data
n=600;
T=5;
p=T;
sigma2=1;
K=512*512;

%Simulate beta1
betai=zeros(n,p);
vaco=iwishrnd(eye(p),14);
save vaco vaco
for i=1:n
    sit=randsample(2,1);
    if sit==1
        betai(i,:)=mvnrnd([1.5 1.5 1 2 2],vaco);
    else
        betai(i,:)=mvnrnd([-1.5 -1.5 -1 -2 -2],vaco);
    end
end
save betai betai

%generate Y and X
for i=1:n
    Yi=zeros(K,T);
    etai=randn(K,1);
    for t=1:T
        posi=zeros(1,T);
        posi(t)=1;
        Xit=kron(ones(K,1),posi);
        Yi(:,t)=etai+(Xit*betai(i,:)')+(1/sqrt(sigma2))*randn(K,1);
        eval(['X_' num2str(i) int2str(t) '=Xit;']);
        eval( ['csvwrite(''X_' num2str(i) int2str(t) '.csv''', ',', 'X_' num2str(i) int2str(t) ')']);
        eval( ['clear X_' num2str(i) int2str(t)]);
        clear Xit
    end
    eval(['Y_' num2str(i) '=Yi;']);
    eval( ['csvwrite(''Y_' num2str(i)  '.csv''', ',', 'Y_' num2str(i) ')']);
    eval( ['clear Y_' num2str(i) ]);
    clear Yi etai
    i
end


indi=[1:K]';
% find to get index & value vectors:
[rind,cind,val] = find(randn(sqrt(K),sqrt(K)));
DD=[rind cind];
clear val rind cind
Cmat=[]; Omega=zeros(K,1);
for k=1:K
    dis=sqrt((DD(k,1)-DD(:,1)).^2+(DD(k,2)-DD(:,2)).^2).^(-2.16);
    dis(isinf(dis))=0;
    Omega(k)=1/sum(dis);
    dis=dis/sum(dis);
    Cmat=[Cmat; [k*ones(sum(dis>=0.01),1) indi(dis>=0.01) dis(dis>=0.01)]];
    k
end
save Cmat Cmat 
save Omega Omega;
clear dis rind cind val DD

%Parameters for priors
asigma=1; bsigma=1;
aalpha=1; balpha=1;
atau=(1/10000); btau=(1/10000);
beta0=randn(1,p); capsigma0=eye(p);
R=20;
npi=40;
pitl=(1/npi)*ones(npi,1);


s=1;
%Read in the data
eval( ['Y_' int2str(s) '= csvread(''Y_' int2str(s) '.csv''' ');']);
eval(['Ys = Y_' num2str(s) ';']);
eval( ['clear Y_' num2str(s) ]);
Xs=[];
for t=1:T
    eval( ['X_' int2str(s) num2str(t) '= csvread(''X_' int2str(s) num2str(t) '.csv''' ');']);
    eval(['Xst = X_' num2str(s) num2str(t) ';']);
    eval( ['clear X_' num2str(s) num2str(t)]);
    Xs=[Xs;Xst];
    clear Xst
end

%Initial values for parameters
asigmatls=1; bsigmatls=1;
aalphatls=1; balphatls=1;
atautls=(1/10000); btautls=(1/10000);
betatls=zeros(p,R);
capsigmatls=zeros(R,p,p);
for r=1:R
    capsigmatls(r,:,:)=eye(p)./10;
end
gama1s=2*ones(1,R-1); gama2s=ones(1,R-1);
wbs=zeros(npi,1);
pitls=pitl;

kappas=rand(1,R);
kappas=kappas/sum(kappas);
xis=mean(Ys,2);
capsis=ones(K,1)./100;

rho=linspace(0,0.9,npi)'; 
load Cmat
so3nr=logdeter(20,20,K,Cmat,rho); 
clear Cmat;
save so3nr so3nr;
clear so3nr
    
[asigmatl,bsigmatl,aalphatl,balphatl,atautl,btautl,betatl,capsigmatl,gama1,gama2,pitl,wb,kappa,xi,capsi]...
        =onlineVB(asigma,bsigma,aalpha,balpha,atau,btau,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,...
        aalphatls,balphatls,atautls,btautls,betatls,capsigmatls,pitls,wbs,kappas,xis,capsis,Ys,Xs,T,K,R,p,rho,s,npi);
asigmatls=asigmatl; bsigmatls=bsigmatl; 
aalphatls=aalphatl; balphatls=balphatl; 
betatls=betatl; capsigmatls=capsigmatl;
gama1s=gama1; gama2s=gama2; 
atautls=atautl; btautls=btautl; wbs=wb;
pitls=pitl;
kappas=kappa;
xis=xi;
capsis=capsi;
clear Ys Xs


s=s+1;
tic;
while s<=n
    %Read in the data
    eval( ['Y_' int2str(s) '= csvread(''Y_' int2str(s) '.csv''' ');']);
    eval(['Ys = Y_' num2str(s) ';']);
    eval( ['clear Y_' num2str(s) ]);
    Xs=[];
    for t=1:T
        eval( ['X_' int2str(s) num2str(t) '= csvread(''X_' int2str(s) num2str(t) '.csv''' ');']);
        eval(['Xst = X_' num2str(s) num2str(t) ';']);
        eval( ['clear X_' num2str(s) num2str(t)]);
        Xs=[Xs;Xst];
        clear Xst
    end

    [asigmatl,bsigmatl,aalphatl,balphatl,atautl,btautl,betatl,capsigmatl,gama1,gama2,pitl,wb,kappa,xi,capsi]...
        =onlineVB(asigma,bsigma,aalpha,balpha,atau,btau,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,...
        aalphatls,balphatls,atautls,btautls,betatls,capsigmatls,pitls,wbs,kappas,xis,capsis,Ys,Xs,T,K,R,p,rho,s,npi);

    asigmatls=asigmatl; bsigmatls=bsigmatl;
    aalphatls=aalphatl; balphatls=balphatl;
    betatls=betatl; capsigmatls=capsigmatl;
    gama1s=gama1; gama2s=gama2;
    atautls=atautl; btautls=btautl; wbs=wb;
    pitls=pitl;
    kappas=kappa;
    xis=xi;
    capsis=capsi;
    s
    s=s+1;
    clear Ys Xs
end
toc;


