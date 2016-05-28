
%use output plus descriptive statistics.mat

burn=1000;
iterused=(burn/thin+1):(round(iter/thin)-1);

%alpha 
alphastat=[mean(precisionout(its,1))' prctile(precisionout(its,1),2.5)' prctile(precisionout(its,1),97.5)']

%sigma
sigma2stat=[mean(1./precisionout(its,2))' prctile(1./precisionout(its,2),2.5)' prctile(1./precisionout(its,2),97.5)']


%True densities
XI=linspace(-4,4,n);
load betai
didit=zeros(p,n);
for l=1:p
    [didi,dodo,widd] = ksdensity(betai(:,l),XI);
    widt(l)=widd;
    didit(l,:)=didi;
end

%Load online VB estimates
load samplvb
didivb=zeros(p,n);
for l=1:p
    [didi,dodo] = ksdensity(samplvb(:,l),XI,'width',widt(l));
    didivb(l,:)=didi;
end

draw=zeros(length(iterused),n,p);
for j=1:length(iterused)
    betastarj=reshape(betastarout(iterused(j),:),R,p);
    weig=vrout(iterused(j),:).*cumprod([1 1-vrout(iterused(j),1:R-1)]);
    nbsampl=randsample(R,n,true,weig);
    nbsampl=histc(nbsampl,1:R);
    sampl=[];
    for r=1:R
        sampl=[sampl;repmat(betastarj(r,:),nbsampl(r),1)];
    end
    for l=1:p
        [didi,dodo] = ksdensity(sampl(:,l),XI,'width',widt(l));
        draw(j,:,l)=didi;
    end
end

drawmean=reshape(mean(draw,1),n,p);
drawp1=reshape(prctile(draw,0.5,1),n,p);
drawp3=reshape(prctile(draw,99.5,1),n,p);
for l=1:p
   subplot(2,3,l)
    plot(XI,didit(l,:),'--b','LineWidth',1.5);
    hold on
    plot(XI,didivb(l,:),'b','LineWidth',1.5);
    hold on
    plot(XI,drawmean(:,l),':k','LineWidth',1.5);
    hold on
    plot(XI,drawp1(:,l),'-.g','LineWidth',3.5);
    hold on
    plot(XI,drawp3(:,l),'-.g','LineWidth',3.5);
    title(['Variable ' , num2str(l)])
    if l==4
        xlim([-5 5])
    else
        xlim([-4 4])
    end
    if l==p
        hleg1 = legend('True','Estimate online VB','Mean estimate MCMC','99% interval estimate MCMC');
        set(hleg1,'Location','Best')
    end
end


