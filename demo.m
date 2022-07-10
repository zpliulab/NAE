% https://github.com/zpliulab/NAE
% zpliu@sdu.edu.cn
clc
clear 
addpath(genpath(pwd))

Plist = [1/3,1/2,2/3,1,2,0];   % Lp

%% data
load 'simu_10.1.mat'

%% rewire
number = 1000;
direct = 1; % 1/0 => Directed/Undirected 
adjlist = rewire(adj,number,1); 

%% NAE
newadjlist = cell(number,1);
newobj=zeros(number,1);
optpvalue = 1;
optobj = 0;
t = [];
Lp = Plist(4);
for i = -8:2
    c.i = i;
    tic
    parfor k = 1:number
        [newgra,obj,par] = refer_NAE(datan,adjlist{k},Lp,c);
        newadjlist{k,1}=newgra;
         newobj(k)=obj(end);
    end
    t(c.i + 9) = toc;
    pvalue = length(find(newobj(1:end)<=newobj(1)))/length(newobj);
    if pvalue < optpvalue
        optpvalue =pvalue;
        optparm = c.i;
    end
    if pvalue == 0.001
        break
    end
end
t = mean (t);

%% p_value
fprintf("\n *********************  Result *************************** \n");
fprintf("   p-value =  %.3f \n",optpvalue);
fprintf("   log_{10}(mu) =  %d \n",optparm);
fprintf("   time(mean)  =  %.4fs \n" ,t);





