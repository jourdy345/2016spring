function [asigmatl,bsigmatl,aalphatl,balphatl,atautl,btautl,betatl,capsigmatl,gama1,gama2,pitl,wb,kappa,xi,capsi]...
        =onlineVB(asigma,bsigma,aalpha,balpha,atau,btau,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,aalphatls,balphatls,atautls,btautls,betatls,capsigmatls,pitls,wbs,kappas,xis,capsis,Ys,Xs,T,K,R,p,rho,s,npi)


%Termination criteria
epsilon = 1e-6;

csp1=1/s;
% csp1=1/(s^0.51);


kappa=kappas;
xi=xis;
capsi=capsis;

%Update subject-specific parameters
prevalue=[kappa(:)' capsi(:)' xi(:)'];
value=zeros(1,length(prevalue));


while mean(abs(prevalue - value)) >= epsilon
    prevalue=[kappa(:)' capsi(:)' xi(:)'];
    
    %kappa
    kappa=zeros(1,R);
    for r=1:R
        sss=0; 
        for t=1:T
            sss=sss-2*Ys(:,t)'*Xs((t-1)*K+1:t*K,:)*betatls(:,r)+2*betatls(:,r)'*Xs((t-1)*K+1:t*K,:)'*xi...
                +betatls(:,r)'*Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*betatls(:,r)+ trace(Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*reshape(capsigmatls(r,:,:),p,p)); 
        end
        if r<R
            kappa(r)=-0.5*(asigmatls./bsigmatls)*(sss)...
                +psi(gama1s(r))-psi(gama1s(r)+gama2s(r))+sum(psi(gama1s(1:(r-1)))-psi(gama1s(1:(r-1))+gama2s(1:(r-1))));
        elseif r==R
            kappa(r)=-0.5*(asigmatls./bsigmatls)*(sss)...
                +sum(psi(gama1s(1:(r-1)))-psi(gama1s(1:(r-1))+gama2s(1:(r-1))));
        end
    end

    kappa=kappa-max(kappa);
    kappa=exp(kappa)/ sum(exp(kappa));

    normalization = false;
    for r=1:R
        if kappa(r) < 1e-10
            kappa(r) = 1e-10; %Avoid domain errors due to limited accuracy of floating point numbers
            normalization = true;
        end
    end
    if normalization == true
        kappa = kappa / sum(kappa); %Normalization
    end
    clear sss ssst sssy sssym
    

    %xi, capsi
    load Omega;
    load Cmat;
    capsi=1./(T*(asigmatls/bsigmatls)+(atautls/btautls)*(1./Omega));
    soom=accumarray(Cmat(:,1),xi(Cmat(:,2)).*Cmat(:,3)); %xi(Cmat(:,1)) replicates elements in xi as many time as they are repeated in Cmat
    soomm=(atautls/btautls)*sum(rho.*pitls)*(1./Omega).*soom;

    s4=zeros(K,1);
    for r=1:R
        for t=1:T
            s4=s4+(asigmatls/bsigmatls)*kappa(r)*(Ys(:,t)-Xs((t-1)*K+1:t*K,:)*betatls(:,r));
        end
    end
    xi=capsi.*(s4+soomm);
    clear s4 xss

    value=[kappa(:)' capsi(:)' xi(:)'];
end
clear value prevalue


%Update subject-invariant parameters
asigmatl=asigmatls; bsigmatl=bsigmatls; aalphatl=aalphatls; balphatl=balphatls; 
betatl=betatls; gama1=gama1s; gama2=gama2s; atautl=atautls; btautl=btautls; pitl=pitls; capsigmatl=capsigmatls; wb=wbs;

prevalue2=[bsigmatl balphatl betatl(:)' gama2(:)' btautl(:)' pitl(:)' capsigmatl(:)'];
value2=zeros(1,length(prevalue2));

sumtx=zeros(p,p); sumyx=zeros(1,p);
for t=1:T
    sumtx=sumtx+Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:);
    sumyx=sumyx+(Ys(:,t)-xi)'*Xs((t-1)*K+1:t*K,:);
end

%asigmatl,
asigmatl=(1-csp1)*asigmatls+csp1*asigma+0.5*T*K;

%aalphatl,
aalphatl=aalpha+R-1;

%gama1tl
gama1=(1-csp1)*gama1s+csp1+kappa(1:(R-1));

%atautl,
atautl=(1-csp1)*atautls+csp1*atau+0.5*K;


while mean(abs(prevalue2 - value2)) >= epsilon
    prevalue2=[bsigmatl balphatl betatl(:)' gama2(:)' btautl(:)' pitl(:)' capsigmatl(:)'];

    %bsigmatl
    ot=0;
    for r=1:R
        sss=0; sssy=0;
        for t=1:T
            sss=sss+trace(Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*reshape(capsigmatl(r,:,:),p,p));
            sssy=sssy+(Ys(:,t)-xi-Xs((t-1)*K+1:t*K,:)*betatl(:,r))'*(Ys(:,t)-xi-Xs((t-1)*K+1:t*K,:)*betatl(:,r));
        end
        ot=ot+kappa(r)*(sssy+sum(capsi)+sss);
    end
    bsigmatl=(1-csp1)*bsigmatls+csp1*bsigma+0.5*ot;
    
    %balphatl
    balphatl=balpha-sum(psi(gama2)-psi(gama1+gama2));
    
    %beta0tl, capsigma0tl
    for r=1:R
        capsigmatlr=inv((1-csp1)*inv(reshape(capsigmatls(r,:,:),p,p))+csp1*inv(capsigma0)+(asigmatl/bsigmatl)*kappa(r)*sumtx);
        capsigmatl(r,:,:)=capsigmatlr;
        betatl(:,r)=(1-csp1)*capsigmatlr*inv(reshape(capsigmatls(r,:,:),p,p))*betatls(:,r)+capsigmatlr*(csp1*beta0*inv(capsigma0)...
            +(asigmatl/bsigmatl)*kappa(r)*sumyx)';
    end
    
    %gama1tl
    gama2=zeros(1,R-1);
    for r=1:(R-1)
        gama2(r)=(1-csp1)*gama2s(r)+(aalphatl/balphatl)*csp1+sum(kappa(r+1:end));
    end

    %btautl
    gamt=xi.*(1./Omega);
    some1=accumarray(Cmat(:,1),gamt(Cmat(:,2)).*Cmat(:,3));
    some=sum(some1.*xi);
    so1=sum(capsi.*(1./Omega))+sum(xi.*(1./Omega).*xi);
    so2=-some;

    btautl=(1-csp1)*btautls+csp1*btau+0.5*so1+0.5*so2*sum(rho.*pitl);

    load so3nr;
    wb=(1-csp1)*wbs-0.5*(atautl/btautl)*so2*rho+0.5*so3nr;
    pitl=wb/(10^(order(wb(end))));
    pitl=exp(pitl)/ sum(exp(pitl));

    normalization = false;
    for ll=1:npi
        if pitl(ll) < 1e-10
            pitl(ll) = 1e-10; %Avoid domain errors due to limited accuracy of floating point numbers
            normalization = true;
        end
    end
    for ll=1:npi
        if normalization == true
            pitl(ll) = pitl(ll) / sum(pitl); %Normalization
        end
    end

    value2=[bsigmatl balphatl betatl(:)' gama2(:)' btautl(:)' pitl(:)' capsigmatl(:)'];
end




