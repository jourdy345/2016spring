function ss=logdeter(p,m,K,Cmat,rho)
%computes log(det(I-rhoC))
vv=zeros(length(rho),p);
for l=1:p
    v=zeros(length(rho),1);
    x=randn(K,1);
    H=x;
    for k=1:m
        H=accumarray(Cmat(:,1),H(Cmat(:,2)).*Cmat(:,3));
        v=(K*(rho.^k).*(x'*H))./k +v;
    end
    v=v./(x'*x);
    vv(:,l)=v;
end
ss=-mean(vv,2);
