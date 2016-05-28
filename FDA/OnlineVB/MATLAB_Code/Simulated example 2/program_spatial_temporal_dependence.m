%Simulate the data
n=400;
T=5;
p=T;
sigma2=1;
K=(512/2)*(512/2);

%Simulate beta1 and Sigmait
betai=zeros(n,p);
vaco=iwishrnd(eye(p),14);
sigmast1=gamrnd(2,3);
sigmast2=gamrnd(4,5);
muT=zeros(n,T);
Sigmait=zeros(n,1);
save vaco vaco

for i=1:n
    sit=randsample(2,1);
    if sit==1
        betai(i,:)=mvnrnd([1.5 1.5 1 2 2],vaco);
        Sigmait(i)=sigmast1;
    else
        betai(i,:)=mvnrnd([-1.5 -1.5 -1 -2 -2],vaco);
        Sigmait(i)=sigmast2;
    end
end
save betai betai
save Sigmait Sigmait

% for l=1:p
%     [didit,dodot] = ksdensity(betai(:,l));
%     subplot(2,3,l)
%     plot(dodot,didit,'--b','LineWidth',1.5);
% end
% 
for i=1:n
    for t=1:T
        if t==1
            muT(i,t)=normrnd(randn,sqrt(1/Sigmait(i)));
        else
            muT(i,t)=normrnd(muT(i,t-1),sqrt(1/Sigmait(i)));
        end
    end
end


%generate Y and X
for i=1:n
    Yi=zeros(K,T);
    etai=randn(K,1);
    for t=1:T
        posi=zeros(1,T);
        posi(t)=1;
        Xit=kron(ones(K,1),posi);
        Yi(:,t)=etai*muT(i,t)+(Xit*betai(i,:)')+(1/sqrt(sigma2))*randn(K,1);
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
Cmat=[]; Omega=zeros(K,1);
for k=1:K
    dis=sqrt((DD(k,1)-DD(:,1)).^2+(DD(k,2)-DD(:,2)).^2).^(-2.16);
    dis(isinf(dis))=0;
    Omega(k)=1/sum(dis);
    dis=dis/sum(dis);
    Cmat=[Cmat; [k*ones(sum(dis>=0.01),1) indi(dis>=0.01) dis(dis>=0.01)]];
end
save Cmat Cmat 
save Omega Omega;
clear dis rind cind val DD




n=400;
T=5;
p=T;
K=(512/2)*(512/2);

%Parameters for priors
asigma=1; bsigma=1;
aalpha=1; balpha=1;
atheta=(1/10000); btheta=(1/10000);
atau=(1/10000); btau=(1/10000);
beta0=randn(1,p); capsigma0=eye(p);
R=20;
npi=40;
pitl=(1/npi)*ones(npi,1);
order=randsample(n,n);
objective=zeros(n,1);


s=1;
%Read in the data
ss=order(s);
eval( ['Y_' int2str(ss) '= csvread(''Y_' int2str(ss) '.csv''' ');']);
eval(['Ys = Y_' num2str(ss) ';']);
eval( ['clear Y_' num2str(ss) ]);
Xs=[];
for t=1:T
    eval( ['X_' int2str(ss) num2str(t) '= csvread(''X_' int2str(ss) num2str(t) '.csv''' ');']);
    eval(['Xst = X_' num2str(ss) num2str(t) ';']);
    eval( ['clear X_' num2str(ss) num2str(t)]);
    Xs=[Xs;Xst];
    clear Xst
end


%Initial values for parameters
asigmatls=1; bsigmatls=1;
aalphatls=1; balphatls=1;
atautls=(1/10000); btautls=(1/10000);
athetatls=(1/10000)*ones(1,R); bthetatls=(1/10000)*ones(1,R);
betatls=zeros(p,R);
capsigmatls=zeros(R,p,p);
for r=1:R
    capsigmatls(r,:,:)=eye(p)./5;
end
gama1s=2*ones(1,R-1); gama2s=ones(1,R-1);
wbs=zeros(npi,1);
pitls=pitl;

kappas=rand(1,R);
kappas=kappas./sum(kappas);
m=3;
xis=zeros(K,m);capsis=zeros(K,m);
for zz=1:m
    xis(:,zz)=mean(Ys,2)/zz;
    capsis(:,zz)=ones(K,1)./100;
end;
lambda1s=randn(m,T);
lambda2s=kron(ones(T,1),eye(m));
lambda01s=zeros(m,1);
lambda02s=eye(m);

rho=linspace(0,0.9,npi)'; 
load Cmat
so3nr=logdeter(20,20,K,Cmat,rho); 
clear Cmat;
save so3nr so3nr;
clear so3nr
    
tic;
[asigmatl,bsigmatl,aalphatl,balphatl,atautl,btautl,athetatl,bthetatl,betatl,capsigmatl,gama1,gama2,pitl,wb,kappa,xi,lambda1,lambda2,lambda01,lambda02,capsi,val]...
        =onlineVB_spatial_temporal_dependence(asigma,bsigma,aalpha,balpha,atau,btau,atheta,btheta,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,aalphatls,balphatls,atautls,btautls,athetatls,bthetatls,betatls,capsigmatls,...
        pitls,wbs,kappas,xis,lambda01s,lambda02s,lambda1s,lambda2s,capsis,Ys,Xs,T,K,R,p,rho,s,npi,n,m);
toc;
asigmatls=asigmatl; bsigmatls=bsigmatl; 
aalphatls=aalphatl; balphatls=balphatl; 
betatls=betatl; capsigmatls=capsigmatl;
gama1s=gama1; gama2s=gama2; 
athetatls=athetatl; bthetatls=bthetatl;
pitls=pitl;
kappas=kappa;
xis=xi;
capsis=capsi;
lambda01s=lambda01;
lambda02s=lambda02;
lambda1s=lambda1;
lambda2s=lambda2;
objective(s)=val;
clear Ys Xs




s=s+1;
tic;
while s<=n
    %Read in the data
    ss=order(s);
    eval( ['Y_' int2str(ss) '= csvread(''Y_' int2str(ss) '.csv''' ');']);
    eval(['Ys = Y_' num2str(ss) ';']);
    eval( ['clear Y_' num2str(ss) ]);
    Xs=[];
    for t=1:T
        eval( ['X_' int2str(ss) num2str(t) '= csvread(''X_' int2str(ss) num2str(t) '.csv''' ');']);
        eval(['Xst = X_' num2str(ss) num2str(t) ';']);
        eval( ['clear X_' num2str(ss) num2str(t)]);
        Xs=[Xs;Xst];
        clear Xst
    end
    
  [asigmatl,bsigmatl,aalphatl,balphatl,atautl,btautl,athetatl,bthetatl,betatl,capsigmatl,gama1,gama2,pitl,wb,kappa,xi,lambda1,lambda2,lambda01,lambda02,capsi,val]...
        =onlineVB_spatial_temporal_dependence(asigma,bsigma,aalpha,balpha,atau,btau,atheta,btheta,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,aalphatls,balphatls,atautls,btautls,athetatls,bthetatls,betatls,capsigmatls,...
        pitls,wbs,kappas,xis,lambda01s,lambda02s,lambda1s,lambda2s,capsis,Ys,Xs,T,K,R,p,rho,s,npi,n,m);


asigmatls=asigmatl; bsigmatls=bsigmatl; 
aalphatls=aalphatl; balphatls=balphatl; 
betatls=betatl; capsigmatls=capsigmatl;
gama1s=gama1; gama2s=gama2; 
athetatls=athetatl; bthetatls=bthetatl;
pitls=pitl;
kappas=kappa;
xis=xi;
capsis=capsi;
lambda01s=lambda01;
lambda02s=lambda02;
lambda1s=lambda1;
lambda2s=lambda2;
objective(s)=val;
s
s=s+1;
    clear Ys Xs
end
toc;



%Discounted Objective function

discoun=zeros(n,1);
for ll=1:n
    disi=1;
    for j=(ll+1):n
        disi=disi*(1-1/(j^0.51));
    end
    discoun(ll)=disi;
end

temp=(psi(asigmatl)-log(bsigmatl))*(asigma-asigmatl)-(asigmatl/bsigmatl)*(bsigma-bsigmatl)+asigma*log(bsigma)-asigmatl*log(bsigmatl)+gammaln(asigmatl)-gammaln(asigma);
valp  = temp;

temp=(psi(aalphatl)-log(balphatl))*(aalpha-aalphatl)-(aalphatl/balphatl)*(balpha-balphatl)+aalpha*log(balpha)-aalphatl*log(balphatl)+gammaln(aalphatl)-gammaln(aalpha);
valp  = valp+ temp;

temp=(psi(atautl)-log(btautl))*(atau-atautl)-(atautl/btautl)*(btau-btautl)+atau*log(btau)-atautl*log(btautl)+gammaln(atautl)-gammaln(atau);
valp  = valp+temp;

temp=(aalphatl/balphatl-1)*(psi(gama2)-psi(gama1+gama2))+psi(aalphatl)-log(balphatl)-(psi(gama1)-psi(gama1+gama2));
valp  =valp  +sum(temp);

temp=(psi(athetatl)-log(bthetatl)).*(atau-athetatl)-(athetatl/bthetatl).*(btau-bthetatl)+atau*log(btau)+gammaln(athetatl)-athetatl.*log(bthetatl)-gammaln(atau);
valp  = valp+sum(temp);

obj=valp+sum(discoun.*objective);
