%Rho
sampl=randsample(length(pitl),1000,true,pitl);
datasamp=zeros(1000,1);
for ii=1:1000
    datasamp(ii)=rho(sampl(ii));
end
destatrho=[mean(datasamp) prctile(datasamp,5) prctile(datasamp,95)]

%alpha
datasamp=gamrnd(aalphatl,1/balphatl,1000,1);
alphastat=[mean(datasamp) prctile(datasamp,2.5) prctile(datasamp,97.5)]

% %sigma
% datasamp=gamrnd(asigmatl,1/bsigmatl,1000,1);
% sigmastat=[mean(datasamp) prctile(datasamp,2.5) prctile(datasamp,97.5)]


%tauj
datasamp=gamrnd(atautl,1/btautl,1000,1);
taustat=[mean(datasamp) prctile(datasamp,2.5) prctile(datasamp,97.5)]

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

plot(weig,'b.--','MarkerSize',15)
xlabel('Component')


nbsampl=randsample(R,(s-1),true,weig);
nbsampl=histc(nbsampl,1:R);
sampl=[];
for r=1:R
    sampl=[sampl;mvnrnd(betatl(:,r)',reshape(capsigmatl(r,:,:),p,p),nbsampl(r))];
end
load betai;

for l=1:p
    [didit,dodot,widd] = ksdensity(betai(:,l));
    [didi,dodo] = ksdensity(sampl(:,l),'width',widd);
    subplot(2,3,l)
    plot(dodo,didi,'b','LineWidth',1.5);
    hold on
    plot(dodot,didit,'--b','LineWidth',1.5);
    title(['Variable ' , num2str(l)])
end


