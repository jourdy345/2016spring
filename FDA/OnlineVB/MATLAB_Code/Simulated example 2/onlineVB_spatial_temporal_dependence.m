function [asigmatl,bsigmatl,aalphatl,balphatl,atautl,btautl,athetatl,bthetatl,betatl,capsigmatl,gama1,gama2,pitl,wb,kappa,xi,lambda1,lambda2,lambda01,lambda02,capsi,val]...
        =onlineVB_spatial_temporal_dependence(asigma,bsigma,aalpha,balpha,atau,btau,atheta,btheta,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,aalphatls,balphatls,atautls,btautls,athetatls,bthetatls,betatls,capsigmatls,...
        pitls,wbs,kappas,xis,lambda01s,lambda02s,lambda1s,lambda2s,capsis,Ys,Xs,T,K,R,p,rho,s,npi,n,m)

%Termination criteria
epsilon = 1e-6;

% csp1=1/s;
csp1=1/(s^0.51);


kappa=kappas;
xi=xis;
capsi=capsis;
lambda01=lambda01s;
lambda02=lambda02s;
lambda1=lambda1s;
lambda2=lambda2s;


%Update subject-specific parameters
prevalue=[kappa(:)' capsi(:)' xi(:)' lambda01(:)' lambda02(:)' lambda1(:)' lambda2(:)'];
value=zeros(1,length(prevalue));


while mean(abs(prevalue - value)) > epsilon
    prevalue=[kappa(:)' capsi(:)' xi(:)' lambda01(:)' lambda02(:)' lambda1(:)' lambda2(:)'];

    %kappa
    kappa=zeros(1,R);
    for r=1:R
        sss=0;saas=0;
        for t=1:T
            sss=sss-2*Ys(:,t)'*Xs((t-1)*K+1:t*K,:)*betatls(:,r)+2*betatls(:,r)'*Xs((t-1)*K+1:t*K,:)'*xi*lambda1(:,t)...
                +betatls(:,r)'*Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*betatls(:,r)+ trace(Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*reshape(capsigmatls(r,:,:),p,p));
            
            if t==1
                saas=saas+(lambda1(:,t)-lambda01)'*(lambda1(:,t)-lambda01)+trace(lambda2(1:m,:)+lambda02);
            else
                saas=saas+(lambda1(:,t)-lambda1(:,t-1))'*(lambda1(:,t)-lambda1(:,t-1))+trace(lambda2((t-1)*m+1:t*m,:)+lambda2((t-2)*m+1:(t-1)*m,:));
            end
        end
        if r<R
            kappa(r)=-0.5*(asigmatls./bsigmatls)*sss -0.5*(athetatls(r)/bthetatls(r))*saas...
                +psi(gama1s(r))-psi(gama1s(r)+gama2s(r))+sum(psi(gama1s(1:(r-1)))-psi(gama1s(1:(r-1))+gama2s(1:(r-1))));
        elseif r==R
            kappa(r)=-0.5*(asigmatls./bsigmatls)*sss -0.5*(athetatls(r)/bthetatls(r))*saas...
                +sum(psi(gama1s(1:(r-1)))-psi(gama1s(1:(r-1))+gama2s(1:(r-1))));
        end
    end
    kappa=kappa-max(kappa-0.01);
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
    did=zeros(m,1);
    for t=1:T
        did=did+diag(lambda2((t-1)*m+1:t*m,:))+lambda1(:,t).*lambda1(:,t);
    end;
    capsi=zeros(K,m);
    for zz=1:m
        capsi(:,zz)=1./((asigmatls./bsigmatls)*did(m)+(atautls/btautls)*(1./Omega));
    end

    s4=zeros(K,m);
    for r=1:R
        for t=1:T
            for ii=1:m
                s4(:,ii)=s4(:,ii)+(asigmatls/bsigmatls)*kappa(r)*lambda1(ii,t)*(Ys(:,t)-Xs((t-1)*K+1:t*K,:)*betatls(:,r));
            end
        end
    end
    if m ==1
       xi=capsi.*s4; 
    elseif m==2
        xi(:,1)=capsi(:,1).*(s4(:,1)-(asigmatls/bsigmatls)*xi(:,2)*sum(lambda1(2,:).*lambda1(1,:)));
        xi(:,2)=capsi(:,2).*(s4(:,2)-(asigmatls/bsigmatls)*xi(:,1)*sum(lambda1(1,:).*lambda1(2,:)));
    else
        xi(:,1)=capsi(:,1).*(s4(:,1)-(asigmatls/bsigmatls)*(xi(:,2)*sum(lambda1(2,:).*lambda1(1,:))+xi(:,3)*sum(lambda1(3,:).*lambda1(1,:))));
        xi(:,2)=capsi(:,2).*(s4(:,2)-(asigmatls/bsigmatls)*(xi(:,1)*sum(lambda1(1,:).*lambda1(2,:))+xi(:,3)*sum(lambda1(3,:).*lambda1(2,:))));
        xi(:,3)=capsi(:,3).*(s4(:,3)-(asigmatls/bsigmatls)*(xi(:,1)*sum(lambda1(1,:).*lambda1(3,:))+xi(:,2)*sum(lambda1(3,:).*lambda1(2,:))));
    end
    clear s4 xss

    %lambdas
    ss1=sum(kappa.*(athetatls./bthetatls));
    ss2=sum(kappa.*(athetatls./bthetatls))*lambda1(:,1);
    lambda02=(1/(ss1+1))*eye(m);
    lambda01=lambda02*ss2;
    
   
    if m==1 
        ded =xi'*xi+sum(capsi);
    elseif m==2
       ded=[xi(:,1)'*xi(:,1)+sum(capsi(:,1)) xi(:,1)'*xi(:,2); xi(:,1)'*xi(:,2) xi(:,2)'*xi(:,2)+sum(capsi(:,2))];
    else
       ded=[xi(:,1)'*xi(:,1)+sum(capsi(:,1)) xi(:,1)'*xi(:,2) xi(:,1)'*xi(:,3); xi(:,1)'*xi(:,2) xi(:,2)'*xi(:,2)+sum(capsi(:,2)) xi(:,2)'*xi(:,3); xi(:,3)'*xi(:,1) xi(:,3)'*xi(:,2) xi(:,3)'*xi(:,3)+sum(capsi(:,3))]; 
    end
    
     for t=1:T
        lambda2((t-1)*m+1:t*m,:)=inv(sum(kappa.*(athetatls./bthetatls))+ded*(asigmatls./bsigmatls));

        suu=zeros(m,1);
        for r=1:R
            suu=suu+kappa(r)*xi'*(Ys(:,t)-Xs((t-1)*K+1:t*K,:)*betatls(:,r));
        end
        if t==1
            lambda1(:,t)=lambda2((t-1)*m+1:t*m,:)*(lambda01.*sum(kappa.*(athetatls./bthetatls))+suu*(asigmatls./bsigmatls));
        else
            lambda1(:,t)=lambda2((t-1)*m+1:t*m,:)*(lambda1(:,t-1)*sum(kappa.*(athetatls./bthetatls))+suu*(asigmatls./bsigmatls));
        end
    end


    value=[kappa(:)' capsi(:)' xi(:)' lambda01(:)' lambda02(:)' lambda1(:)' lambda2(:)'];

   %mean(abs(prevalue - value)) 
end
clear value prevalue


%Update subject-invariant parameters
asigmatl=asigmatls; bsigmatl=bsigmatls; aalphatl=aalphatls; balphatl=balphatls; athetatl=athetatls; bthetatl=bthetatls;
betatl=betatls; gama1=gama1s; gama2=gama2s; atautl=atautls; btautl=btautls; pitl=pitls; capsigmatl=capsigmatls; wb=wbs; 

prevalue2=[bsigmatl balphatl betatl(:)' gama2(:)' btautl(:)' bthetatl(:)' pitl(:)' capsigmatl(:)'];
value2=zeros(1,length(prevalue2));

sumtx=zeros(p,p); sumyx=zeros(1,p);
for t=1:T
    sumtx=sumtx+Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:);
    sumyx=sumyx+(Ys(:,t)-xi*lambda1(:,t))'*Xs((t-1)*K+1:t*K,:);
end

%asigmatl,
asigmatl=(1-csp1)*asigmatls+csp1*asigma+0.5*T*K;

%aalphatl,
aalphatl=aalpha+R-1;

%gama1tl
gama1=(1-csp1).*gama1s+csp1+kappa(:,1:(R-1));

%atautl,
atautl=(1-csp1)*atautls+csp1*atau+0.5*K;

%athetatl
athetatl=(1-csp1)*athetatls+csp1*atheta+0.5*T*kappa(1,:);


while mean(abs(prevalue2 - value2)) >= epsilon
    prevalue2=[bsigmatl balphatl betatl(:)' gama2(:)' btautl(:)' bthetatl(:)' pitl(:)' capsigmatl(:)'];

    %bsigmatl
    ot=0;
    if m==1 
        dod =xi'*xi+sum(capsi);
    elseif m==2
       dod=[xi(:,1)'*xi(:,1)+sum(capsi(:,1)) xi(:,1)'*xi(:,2); xi(:,1)'*xi(:,2) xi(:,2)'*xi(:,2)+sum(capsi(:,2))];
    else
       dod=[xi(:,1)'*xi(:,1)+sum(capsi(:,1)) xi(:,1)'*xi(:,2) xi(:,1)'*xi(:,3); xi(:,1)'*xi(:,2) xi(:,2)'*xi(:,2)+sum(capsi(:,2)) xi(:,2)'*xi(:,3); xi(:,3)'*xi(:,1) xi(:,3)'*xi(:,2) xi(:,3)'*xi(:,3)+sum(capsi(:,3))]; 
    end
    for r=1:R
        sss=0; sssy=0;ssso=0;
        for t=1:T
            sss=sss+trace(Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*reshape(capsigmatl(r,:,:),p,p));
            sssy=sssy+sum((Ys(:,t)-xi*lambda1(:,t)-Xs((t-1)*K+1:t*K,:)*betatl(:,r)).^2);
            ssso=ssso+trace(dod*lambda2((t-1)*m+1:t*m,:))+lambda1(:,t)'*dod*lambda1(:,t)-lambda1(:,t)'*(xi'*xi)*lambda1(:,t);
        end
        ot=ot+kappa(r)*(sssy+ssso+sss);
    end
    bsigmatl=(1-csp1)*bsigmatls+csp1*bsigma+0.5*ot;

    %balphatl
    balphatl=balpha-sum(psi(gama2)-psi(gama1+gama2));

    %bthetatl
    bthetatl=zeros(1,R);
    for r=1:R
        su=0;
        for t=1:T
            if t==1
                su=su+kappa(r)*((lambda1(:,t)-lambda01)'*(lambda1(:,t)-lambda01)+trace(lambda2(1:m,:)+lambda02));
            else
                su=su+kappa(r)*(lambda1(:,t)-lambda1(:,t-1))'*(lambda1(:,t)-lambda1(:,t-1))+trace(lambda2((t-1)*m+1:t*m,:)+lambda2((t-2)*m+1:(t-1)*m,:));
            end
        end
        bthetatl(r)=(1-csp1)*bthetatls(r)+csp1*btheta+0.5*su;
    end


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
    ssp=0;some=0;
    for hh=1:m
        gamt=xi(:,hh).*(1./Omega);
        some1=accumarray(Cmat(:,1),gamt(Cmat(:,2)).*Cmat(:,3));
        some=some+sum(some1.*xi(:,hh));
        ssp=ssp+sum(xi(:,hh).*(1./Omega).*xi(:,hh));
    end
    if m==1
        so1=m*sum(capsi.*(1./Omega))+ssp;
    elseif m==2
        so1=sum(capsi(:,1).*(1./Omega)+capsi(:,2).*(1./Omega))+ssp;
    else
        so1=sum(capsi(:,1).*(1./Omega)+capsi(:,2).*(1./Omega)+capsi(:,3).*(1./Omega))+ssp;
    end;
    so2=-some;

    btautl=(1-csp1)*btautls+csp1*btau+0.5*so1+0.5*so2*sum(rho.*pitl);
    
    
    %pitl
    %previously,
    load so3nr;
    wb=(1-csp1)*wbs-0.5*(atautl/btautl)*so2*rho+0.5*so3nr;
    pitl=wb-max(wb-0.01);
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

    value2=[bsigmatl balphatl betatl(:)' gama2(:)' btautl(:)' bthetatl(:)' pitl(:)' capsigmatl(:)'];
    %mean(abs(prevalue2 - value2))
end


%Calculate objective function

% discoun=1;
% for l=(s+1):n
%     discoun=discoun*(1-1/(l^0.51));
% end
    

% temp=(psi(asigmatl)-log(bsigmatl))*(asigma-asigmatl)-(asigmatl/bsigmatl)*(bsigma-bsigmatl)+asigma*log(bsigma)-asigmatl*log(bsigmatl)+gammaln(asigmatl)-gammaln(asigma);
% val  = temp;
% 
% temp=(psi(aalphatl)-log(balphatl))*(aalpha-aalphatl)-(aalphatl/balphatl)*(balpha-balphatl)+aalpha*log(balpha)-aalphatl*log(balphatl)+gammaln(aalphatl)-gammaln(aalpha);
% val  = val+ temp;
% 
% temp=(psi(atautl)-log(btautl))*(atau-atautl)-(atautl/btautl)*(btau-btautl)+atau*log(btau)-atautl*log(btautl)+gammaln(atautl)-gammaln(atau);
% val  = val+temp;
% 
% temp=(aalphatl/balphatl-1)*(psi(gama2)-psi(gama1+gama2))+psi(aalphatl)-log(balphatl)-(psi(gama1)-psi(gama1+gama2));
% val  =val  +sum(temp);
% 
% temp=(psi(athetatl)-log(bthetatl)).*(atau-athetatl)-(athetatl/bthetatl).*(btau-bthetatl)+atau*log(btau)+gammaln(athetatl)-athetatl.*log(bthetatl)-gammaln(atau);
% val  = val+sum(temp);

temp=0;
for r=1:(R-1)
    temp=temp+kappa(r)*((psi(gama1(r))-psi(gama1(r)+gama2(r)))+sum(psi(gama2(1:(r-1)))-psi(gama1(1:(r-1))+gama2(1:(r-1)))))...
        -kappa(r)*log(kappa(r));
end

temp=temp+kappa(R)*(sum(psi(gama2(1:(R-1)))-psi(gama1(1:(R-1))+gama2(1:(R-1)))))...
    -kappa(R)*log(kappa(R));
val  = temp;


temp=0;
for t=1:T
    if t==1
        temp=temp-0.5*((lambda1(:,t)-lambda01)'*(lambda1(:,t)-lambda01)+trace(lambda2(1:m,:)+lambda02))*sum(kappa.*(athetatl./bthetatl))...
            +0.5*sum(kappa.*(psi(athetatl)-log(bthetatl)-log(2*pi)))-0.5*log(det(lambda2((t-1)*m+1:t*m,:)));
    else
        temp=temp-0.5*(lambda1(:,t)-lambda1(:,t-1))'*(lambda1(:,t)-lambda1(:,t-1))+trace(lambda2((t-1)*m+1:t*m,:)+lambda2((t-2)*m+1:(t-1)*m,:))*sum(kappa.*(athetatl./bthetatl))...
            +0.5*sum(kappa.*(psi(athetatl)-log(bthetatl)-log(2*pi)))-0.5*log(det(lambda2((t-1)*m+1:t*m,:)));
    end
end
val  =val + temp;

if m==1
    temp=-0.5*sum(log(capsi))+0.5*m*K*(psi(atautl)-log(btautl))-0.5*sum(pitl.*so3nr)-0.5*(atautl/btautl)*(so1+so2*sum(rho.*pitl));
    val=val+temp;
elseif m==2
    temp=-0.5*sum(log(capsi(:,1)))-0.5*sum(log(capsi(:,2)))+0.5*m*K*(psi(atautl)-log(btautl))-0.5*sum(pitl.*so3nr)-0.5*(atautl/btautl)*(so1+so2*sum(rho.*pitl));
    val=val+temp;
else
    temp=-0.5*sum(log(capsi(:,1)))-0.5*sum(log(capsi(:,2)))-0.5*sum(log(capsi(:,3)))+0.5*m*K*(psi(atautl)-log(btautl))-0.5*sum(pitl.*so3nr)-0.5*(atautl/btautl)*(so1+so2*sum(rho.*pitl));
    val=val+temp;
end

ot=0;
    if m==1 
        dod =xi'*xi+sum(capsi);
    elseif m==2
       dod=[xi(:,1)'*xi(:,1)+sum(capsi(:,1)) xi(:,1)'*xi(:,2); xi(:,1)'*xi(:,2) xi(:,2)'*xi(:,2)+sum(capsi(:,2))];
    else
       dod=[xi(:,1)'*xi(:,1)+sum(capsi(:,1)) xi(:,1)'*xi(:,2) xi(:,1)'*xi(:,3); xi(:,1)'*xi(:,2) xi(:,2)'*xi(:,2)+sum(capsi(:,2)) xi(:,2)'*xi(:,3); xi(:,3)'*xi(:,1) xi(:,3)'*xi(:,2) xi(:,3)'*xi(:,3)+sum(capsi(:,3))]; 
    end
for r=1:R
    sss=0; sssy=0;ssso=0;
    for t=1:T
        sss=sss+trace(Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*reshape(capsigmatl(r,:,:),p,p));
        sssy=sssy+(Ys(:,t)-xi*lambda1(:,t)-Xs((t-1)*K+1:t*K,:)*betatl(:,r))'*(Ys(:,t)-xi*lambda1(:,t)-Xs((t-1)*K+1:t*K,:)*betatl(:,r));
        ssso=ssso+trace(dod*lambda2((t-1)*m+1:t*m,:))+lambda1(:,t)'*dod*lambda1(:,t)-lambda1(:,t)'*(xi'*xi)*lambda1(:,t);
    end
    ot=ot+kappa(r)*(sssy+sssy+ssso+sss);
end
val  =val + (-0.5*(asigmatl/bsigmatl)*ot+0.5*T*(psi(asigmatl)-log(bsigmatl)-log(2*pi)));





