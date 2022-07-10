function [result,obj] = L2norm(Xr,Xs,xx,lambda)
genenum = sqrt(size(xx,1));
II = [eye(size(xx,1)) -eye(size(xx,1))];
Dmat = II' * (Xr' * Xr + lambda) *II ;
dvec = ((Xs' * Xr * II))';

  Ieq = find(xx==0);
  if(isempty(Ieq))
    E = [];
    F = [];
  else
    E = zeros(length(Ieq),size(xx,1));
    j = 0;
    
    for i = 1:size(xx,1)
      if(xx(i,1)==0)
        j = j+1;
        E(j,i) = 1;
      end
    end
    F = zeros(size(E,1),1);
  end
 Ineg = find(xx==-1);
  if(isempty(Ineg))
    N = [];
    M = [];
  else
    N = zeros(length(Ineg),size(xx,1));
    j = 0;
    for i = 1:size(xx,1)
      if xx(i,1) == -1
        j = j+1;
        N(j,i) = -1;
      end
    end
    M = zeros(size(N,1),1);
  end
  if ~isempty(N)
    N = N * II;
  end
  if ~isempty(E)
    E = E * II;
  end
  Amat  = [E; N; eye(2*size(xx,1))];
  G = zeros(2*size(xx,1),1);
  bvec  = [F;M;G];
opts = optimset('Display','off');
[x,obj] = quadprog(Dmat,-1*dvec,-1*Amat(length(Ieq)+1:end,:),bvec(length(Ieq)+1:end,:),...
                    Amat(1:length(Ieq),:),bvec(1:length(Ieq),:),[],[],[],opts);
result = x(1:size(x,1)/2) - x(size(x,1)/2+1:end);
result = reshape(result,genenum,genenum)';
result(abs(result)<1e-6)=0;
obj = abs(obj);
end