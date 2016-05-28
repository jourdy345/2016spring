function [asigmatl,bsigmatl,aalphatl,balphatl,betatl,capsigmatl,gama1,gama2,kappa]...
        =onlineVB_noeta(asigma,bsigma,aalpha,balpha,beta0,capsigma0,gama1s,gama2s,asigmatls,bsigmatls,aalphatls,balphatls,betatls,capsigmatls,kappas,Ys,Xs,T,K,R,p,s)


%Termination criteria
epsilon = 1e-5;

csp1=1/s;
% csp1=1/(s^0.51);
% csp1=0;



kappa=kappas;

%Update subject-specific parameters
prevalue=[kappa(:)'];
value=zeros(1,length(prevalue));


while mean(abs(prevalue - value)) >= epsilon
    prevalue=[kappa(:)'];
    
    %kappa
    kappa=zeros(1,R);
    for r=1:R
        sss=0; 
        for t=1:T
            sss=sss-2*Ys(:,t)'*Xs((t-1)*K+1:t*K,:)*betatls(:,r)+...
                +betatls(:,r)'*Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*betatls(:,r); 
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
    

    value=[kappa(:)'];
end
clear value prevalue


%Update subject-invariant parameters
asigmatl=asigmatls; bsigmatl=bsigmatls; aalphatl=aalphatls; balphatl=balphatls; 
betatl=betatls; gama1=gama1s; gama2=gama2s; capsigmatl=capsigmatls; 

prevalue2=[bsigmatl balphatl betatl(:)' gama2(:)' capsigmatl(:)'];
value2=zeros(1,length(prevalue2));

sumtx=zeros(p,p); sumyx=zeros(1,p);
for t=1:T
    sumtx=sumtx+Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:);
    sumyx=sumyx+(Ys(:,t))'*Xs((t-1)*K+1:t*K,:);
end

%asigmatl,
asigmatl=(1-csp1)*asigmatls+csp1*asigma+0.5*T*K;

%aalphatl,
aalphatl=aalpha+R-1;

%gama1tl
gama1=(1-csp1)*gama1s+csp1+kappa(1:(R-1));


while mean(abs(prevalue2 - value2)) >= epsilon
    prevalue2=[bsigmatl balphatl betatl(:)' gama2(:)' capsigmatl(:)'];

    %bsigmatl
    ot=0;
    for r=1:R
        sss=0; sssy=0;
        for t=1:T
            sss=sss+trace(Xs((t-1)*K+1:t*K,:)'*Xs((t-1)*K+1:t*K,:)*reshape(capsigmatl(r,:,:),p,p));
            sssy=sssy+(Ys(:,t)-Xs((t-1)*K+1:t*K,:)*betatl(:,r))'*(Ys(:,t)-Xs((t-1)*K+1:t*K,:)*betatl(:,r));
        end
        ot=ot+kappa(r)*(sssy+sss);
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


    value2=[bsigmatl balphatl betatl(:)' gama2(:)' capsigmatl(:)'];
end




