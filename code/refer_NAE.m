function [result,obj,optpar] = refer_NAE(datdata,adj,p,c)
exp_y = datdata(:,2:end);
exp_x = datdata(:,1:end-1);
Xr = exp_x';
Xs = exp_y';
xx = adj';
lambdamax = 2*norm((Xr' * Xs),Inf);
lambda = lambdamax*10.^(-5:0);
optloss = 1e6;
optpar=[];
%%
if p < 2 % P <2 
    for opt = 1:length(lambda)  
        for i = c.i
            opts.mu1 = 10^(i);
            opts.nk = 1;
            opts.lambda = lambda(opt);
            opts.p = p;
           [~,loss] = NAE(Xr,Xs,xx,opts);
           kk = 1;
           while (loss(end) >1e3) && (kk <= 3)
                opts.nk = 0.1*opts.nk;
                [~,loss] = NAE(Xr,Xs,xx,opts);
                kk = kk+1;
           end
            if loss(end) < optloss
                optloss = loss(end);
                optpar.mu1 = 10^(i);
                optpar.nk = opts.nk;
                optpar.lambda = lambda(opt);
                optpar.p = p;
            end
        end
    end
    if isfield(optpar,'mu1')
        [result,obj] = NAE(Xr,Xs,xx,optpar);   
    else
        obj = Inf;
        result=[];
    end
else % P = 2 
      
    for opt = 1:length(lambda)  
        opts.lambda = lambda(opt);
        opts.p = p;
        Xr = kron(eye(length(adj)),exp_x)';
        Xs = reshape(exp_y',[],1);
        xx = reshape(adj',length(adj)^2,1);
       [~,loss] = L2norm(Xr,Xs,xx,lambda(opt));
        if loss(end) < optloss
            optloss = loss(end);
            optpar.lambda = lambda(opt);
            optpar.p = p;
        end
    end
     if isfield(optpar,'lambda')
        [result,obj] = L2norm(Xr,Xs,xx,lambda(opt));

    else
        obj = Inf;
        result=[];
    end
    
end
end


