function [result,obj,optpar] = refer_NAE(dat,adj,p,c)
datdata =dat.data;
exp_y = datdata(:,2:end);
exp_x = datdata(:,1:end-1);
Xr = exp_x';
Xs = exp_y';
xx = adj';
lambdamax = 2*norm((Xr' * Xs),Inf);
lambda = lambdamax*10.^(-4:0);
optloss = 1e6;
optpar=[];
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
end


