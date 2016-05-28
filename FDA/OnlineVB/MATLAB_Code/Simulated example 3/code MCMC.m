n=400;
T=5;
p=T;
sigma2=1;
K=(50)*(50);

%Load the same data used for online VB
Y=zeros(n,K,T); X=zeros(T,K,T);
for ss=1:n
    eval( ['Y_' int2str(ss) '= csvread(''Y_' int2str(ss) '.csv''' ');']);
    eval(['Ys = Y_' num2str(ss) ';']);
    eval( ['clear Y_' num2str(ss) ]);
    Y(ss,:,:)=Ys;
end

for t=1:T
    eval( ['X_' int2str(ss) num2str(t) '= csvread(''X_' int2str(ss) num2str(t) '.csv''' ');']);
    eval(['Xst = X_' num2str(ss) num2str(t) ';']);
    eval( ['clear X_' num2str(ss) num2str(t)]);
    X(t,:,:)=Xst;
    clear Xst
end
    

%Parameters for priors
asigma=1; bsigma=1;
aalpha=1; balpha=1;
ap=1; bp=1;
atau=(1/10000); btau=(1/10000);
beta00=randn(1,p); capsigma00=eye(p)/10;
R=20;

nrun = 5000;  %number runs
burn = 1000;  %burn-in
thin = 5;  %thinning


% --- Initial values --- %
alpha0 = 1;                        % alpha
zi = unidrnd(R,n,1);     % component indicator (1,2,...,M)
betastar = randn(R,p);        % atoms shared ny all mixture distributions
sigma20=1;
beta0=randn(n,p);
ps=1;

% --- Define output files --- %
precisionout = zeros((nrun)/thin,2);            % [alpha, sigma20]
vrout=zeros((nrun)/thin,R);
betastarout=zeros((nrun)/thin,R*p);



iter=1;

% -- Gibbs sampling -- %
while iter<=nrun

    %     tic;
    % -- update vr -- %
    %     tic;
    for r=1:R-1
        vr0(r)=betarnd(1+sum(zi==r),alpha0+sum(zi>r));
    end
    vr0(R)=1;
    %     fprintf('vr: Elapsed Time = %7.2f Seconds\n',toc)

    % -- update zi -- %
    %     tic;
    pzi=zeros(n,R);
    pitk=vr0.*cumprod([1 1-vr0(1:R-1)]);
    for r=1:R
        lik=zeros(n,K,T);
        for t=1:T
            lik(:,:,t)=log(normpdf(Y(:,:,t),betastar(r,t),sqrt(1/ps)));
        end
        pzi(:,r)=log(repmat(pitk(r),1,n))+(sum(reshape(sum(lik,2),n,T),2)');
    end
    pzi=pzi-repmat(max(pzi,[],2),1,R);
    pzi=exp(pzi)./ repmat(sum(exp(pzi),2),1,R);
    zi=mymnrnd(1,pzi).*repmat(1:R,n,1);
    zi=sum(zi,2);
    %     fprintf('zi: Elapsed Time = %7.2f Seconds\n',toc)


    % -- Update unique values betastar -- %
    %     tic;
    for r =1:R
        sss=zeros(p,p); ssmu=zeros(sum(zi==r),p);
        for t=1:T
            Xtre=reshape(X(t,:,:),K,p);
            sss=sss+Xtre'*Xtre;
            ssmu=ssmu+ps*Y(zi==r,:,t)*Xtre;
        end
        sigr=inv(inv(capsigma00)+ps*sum(zi==r)*sss);
        mur=(sum(ssmu,1))*sigr;
        betastarr=mvnrnd(mur,sigr);
        betastar(r,:)=betastarr;
        beta0(zi==r,:)=repmat(betastarr,sum(zi==r),1);
    end
    clear sigk  mu1k deltastark Xtre
    %     fprintf('Deltastar: Elapsed Time = %7.2f Seconds\n',toc)


    % -- update 1/sigma2, ps
    %     tic;
    mupsa=0;
    for t=1:T
        mupsa= mupsa+sum(sum((Y(:,:,t)'-reshape(X(t,:,:),K,p)*beta0').^2));
    end
    ps =gamrnd(ap+n*T*K/2,1/(bp+0.5*mupsa));
    %     fprintf('ps: Elapsed Time = %7.2f Seconds\n',toc)


    % -- Update alpha -- %
    %     tic;
    xx=betarnd(alpha0+1,n);
    phii=(aalpha+length(unique(zi))-1)/(n*(balpha-log(xx))+aalpha+length(unique(zi)-1));
    rr=simdiscrete([phii 1-phii],1);
    if rr==1
        alpha0=gamrnd(aalpha+ length(unique(zi)),1./(balpha-sum(log(xx))));
    else
        alpha0=gamrnd(aalpha+ length(unique(zi))-1,1./(balpha-sum(log(xx))));
    end
    clear xx rr phii
    %     fprintf('alpha: Elapsed Time = %7.2f Seconds\n',toc)


    % -- save sampled values (after thinning) -- %
    if mod(iter,thin)==0
        precisionout(iter/thin,:)=[alpha0 ps];
        vrout(iter/thin,:)=vr0(:)';
        betastarout(iter/thin,:)=betastar(:)';
    end
    [iter alpha0 ps]
    iter=iter+1;
    %         toc;
end




