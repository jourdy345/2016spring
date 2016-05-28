function r=lmvnpdf(a,b,C)
        [nr,nc]=size(a);
%         r=-(nr/2)*log((2*pi))-0.5*log(det(C))-0.5*sum(sum(((a-b)*(inv(C))*(a-b)')));
        r=-0.5*sum(sum(((a-b)*(inv(C))*(a-b)')));
    end