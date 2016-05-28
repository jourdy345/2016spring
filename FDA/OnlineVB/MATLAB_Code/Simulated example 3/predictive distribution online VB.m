%mixing distribution DP
% R=20;
weig=zeros(1,R);
for r=1:(R-1)
    if r==1
        weig(r)=(gama1(r)/(gama1(r)+gama2(r)));
    else
        weig(r)=(gama1(r)/(gama1(r)+gama2(r)))*prod(gama2(1:(r-1))./(gama1(1:(r-1))+gama2(1:r-1)));
    end
end
weig(R)=1-sum(weig(1:(R-1)));

% plot(weig,'b.--','MarkerSize',15)
% xlabel('Component')


%graph betai with width

nbsampl=randsample(R,n,true,weig);
nbsampl=histc(nbsampl,1:R);
sampl=[];
for r=1:R
    sampl=[sampl;mvnrnd(betatl(:,r)',reshape(capsigmatl(r,:,:),p,p),nbsampl(r))];
end
load betai
for l=1:p
    [didit,dodot,widd] = ksdensity(betai(:,l));
    [didi,dodo] = ksdensity(sampl(:,l),'width',widd);
    subplot(2,3,l)
    plot(dodo,didi,'b','LineWidth',1.5);
    hold on
    plot(dodot,didit,'--b','LineWidth',1.5);
    title(['Variable ' , num2str(l)])
    if l==p
        hleg1 = legend('Estimate online VB','True');
        set(hleg1,'Location','Best')
    end
end

%Save online VB estimates to compare with those obtained using MCMC
save samplvb sampl


