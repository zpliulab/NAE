function adjlist = rewire(adj,number,method)
%% rewire    
% number: Number of rewiring
% method: = 1,   directed network;
%         other, undirected network.
if nargin  < 2
    number = 1000;
end
if nargin  < 3
    method = 1;
end
adjlist =cell(1,1);
adjlist{1,1}=adj;
k = 1;
new_adj = adj;
time = 1;
rewire = 2;
while rewire
    if rewire == 2
        eq = 0;
        if method == 1
            new_adj=dir_generate_srand(adj);
        else
            new_adj=sym_generate_srand(adj);
        end
        for i = 1 : length(adjlist)
            if isequal(new_adj, adjlist{i,1}) 
                eq = 1;
                time = time+1;
                break
            end
        end
        if eq == 0
            k = k+1;
            new_adj(new_adj>1) = 1;
            adjlist{k,1} = new_adj;
            time = 1;
        end
        if k >= number
           rewire = 0; 
        end
        if time >=1e5
            rewire = 1;
        end
    else
       if method == 1
            new_adj=dir_generate_srand(adj);
        else
            new_adj=sym_generate_srand(adj);
        end
        if ~isequal(new_adj,adj) 
             k = k+1;
             new_adj(new_adj>1) = 1;
            adjlist{k,1} = new_adj;            
        end
        if k >= number
           rewire = 0; 
        end
    end   
end
end