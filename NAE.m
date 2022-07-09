function [result,obj] = NAE(Xr,Xs,xx,opts)
mu1 = opts.mu1;
nk = opts.nk;
p = opts.p;
lambda = opts.lambda;
display = 0;
%================================  不等式条件构建 =============================== 
  Et = zeros(size(xx));
  Et ( xx == 0) = 1;
  St = zeros(size(xx));
  St ( xx == -1) = 1;
%================================   
Xrt =  Xr;
xxt = xx;
X = Xrt;
B = Xs;
A = xxt;
%% ================================  
Q =  A;
I  = eye(size(A));
epl = zeros(size(A))+1e-5;
maxIter = 500; 
C1 = sparse(zeros(size(A)));
C2 = sparse(zeros(size(A)));
C3 = sparse(zeros(size(A)));
XB = 2*X'*B;
XTX = sparse(X'*X);
for ii=1:maxIter
    %================================  update A =============================== 
    A = A - nk*(2*XTX*A - XB + mu1*((St.*I)*(St.*A+epl+C1/mu1)+(Et.*I)*(Et.*A+C2/mu1) + A -Q) + C1 + C2 + C3);
    %================================  update Q ===============================
    c = A + C1 / mu1;
    if p == 1
        Q = prox_1(c,lambda/mu1);
    elseif p==0
        Q = prox_0(c,lambda/mu1);   
    else
        Q = prox_p(c,lambda/mu1,p);        
    end
    %================================ update C ===============================
    C1 = max(C1 + mu1 * (St.*A+epl),0);
    C2 = C2 + mu1 * (Et.*A);
    C3 = C3 + mu1 * (A-Q);
    mu1 = mu1 * 1.1;
    SA1 = St.*A;
    C1(SA1~=0) = 0;
    %================================ update sita ===============================   
    epl = max(-(St.*A + C1/ mu1) ,0);

    %================================ obj ===============================
    obj(ii) = norm(X*A-B,2)^2;
    if display == 1
        fprintf('itr: %4d\tfval: %epl\tfeasi:%.1e\n', ii, f,obj(ii));
    end
     if ii >50 && (obj(ii)>obj(1))
        break
     end
    
     if ii >1
         if abs(obj(ii)-obj(ii-1))<1e-5
            break
         end
         if (obj(ii)-obj(ii-1))>1e3
            break
         end
     end
end
result = A;
result(abs(result)<1e-3)=0;
end

function Q = prox_0(rho, ccc)   
    Q = rho;
    Q(rho.^2<(2/(1/ccc))) = 0;
end

function Q = prox_p(D, alpha,p)   
    DD = max(abs(D) - alpha*(abs(D).^(p-1)),0);
    Q = DD.*D./abs(D);
end

function y = prox_1(x, mu)
y = max(abs(x) - mu, 0);
y = sign(x) .* y;
end