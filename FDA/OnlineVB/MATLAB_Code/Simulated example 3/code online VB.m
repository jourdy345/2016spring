%Simulate the data
n=400;
T=5;
p=T;
sigma2=1;
K=(50)*(50);

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

% for l=1:p
%     [didit,dodot] = ksdensity(betai(:,l));
%     subplot(2,3,l)
%     plot(dodot,didit,'--b','LineWidth',1.5);
% end
% 
%generate Y and X
for i=1:n
    Yi=zeros(K,T);
    for t=1:T
        posi=zeros(1,T);
        posi(t)=1;
        Xit=kron(ones(K,1),posi);
        Yi(:,t)=(Xit*betai(i,:)')+(1/sqrt(sigma2))*randn(K,1);
        eval(['X_' num2str(i) int2str(t) '=Xit;']);
        eval( ['csvwrite(''X_' num2str(i) int2str(t) '.csv''', ',', 'X_' num2str(i) int2str(t) ')']);
        eval( ['clear X_' num2str(i) int2str(t)]);
        clear Xit
    end
    eval(['Y_' num2str(i) '=Yi;']);
    eval( ['csvwrite(''Y_' num2str(i)  '.csv''', ',', 'Y_' num2str(i) ')']);
    eval( ['clear Y_' num2str(i) ]);
    clear Yi 
    i
end


%Parameters for priors
asigma=1; bsigma=1;
aalpha=1; balpha=1;
beta0=randn(1,p); capsigma0=eye(p);
R=20;
order=randsample(n,n);

tic;
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
betatls=zeros(p,R);
capsigmatls=zeros(R,p,p);
for r=1:R
    capsigmatls(r,:,:)=eye(p)./10;
end
gama1s=2*ones(1,R-1); gama2s=ones(1,R-1);

kappas=rand(1,R);
kappas=kappas/sum(kappas);

    
% tic;
[asigmatl,bsigmatl,aalphatl,balphatl,betatl,capsigmatl,gama1,gama2,kappa]...
        =onlineVB_noeta(asigma,bsigma,aalpha,balpha,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,...
        aalphatls,balphatls,betatls,capsigmatls,kappas,Ys,Xs,T,K,R,p,s);
% toc;
asigmatls=asigmatl; bsigmatls=bsigmatl; 
aalphatls=aalphatl; balphatls=balphatl; 
betatls=betatl; capsigmatls=capsigmatl;
gama1s=gama1; gama2s=gama2; 
kappas=kappa;
clear Ys Xs




s=s+1;
% tic;
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
    
    %Initial values
    

    %     tic;
[asigmatl,bsigmatl,aalphatl,balphatl,betatl,capsigmatl,gama1,gama2,kappa]...
        =onlineVB_noeta(asigma,bsigma,aalpha,balpha,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,...
        aalphatls,balphatls,betatls,capsigmatls,kappas,Ys,Xs,T,K,R,p,s);
    %     toc;

    asigmatls=asigmatl; bsigmatls=bsigmatl;
    aalphatls=aalphatl; balphatls=balphatl;
    betatls=betatl; capsigmatls=capsigmatl;
    gama1s=gama1; gama2s=gama2;
    kappas=kappa;
    s
    s=s+1;
    clear Ys Xs
end
toc;


