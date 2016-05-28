function sample=simdiscrete(p,m)
%General discrete distribution  
%Suppose probabilities p=[p(1) ... p(n)] for the values [1:n] are given, sum(p)=1 and the components p(j) 
%are nonnegative. To generate a random sample of size m from this distribution imagine that the interval 
%(0,1) is divided into intervals with the lengths p(1),...,p(n). Generate a uniform number rand, if this 
%number falls in the jth interval give the discrete distribution the value j. Repeat m times. 
% if sum(p)~=1.0000, disp('Error: the sum of the probabilities must be 1'), return, end;
[r,n]=size(p);
uni=rand(1,m);
cumprob=[0 cumsum(p)];
sample=zeros(1,m);
for j=1:n
  ind=find((uni>cumprob(j)) & (uni<=cumprob(j+1)));
  sample(ind)=j;
end





