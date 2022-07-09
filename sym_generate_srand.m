function srand=sym_generate_srand(s1,ntry)
% Syntax: 
% srand=sym_generate_srand(s1)
% s1 - the adjacency matrix of an undirected network  
% ntry - (optional) the number of rewiring steps. If none is given ntry=4*(# of edges in the network)
% Output: srand - the adjacency matrix of a randomized network with the same set of in- and out-degrees as the original one 


nrew=0;
srand=s1;
[i_srand,j_srand]=find(srand);
aux=find(i_srand>j_srand);
i_srand=i_srand(aux);
j_srand=j_srand(aux);
Ne=length(i_srand);
if (nargin < 2) 
    ntry = 4*Ne; 
end

for i=1:ntry      
   e1 = 1+floor(Ne*rand);
   e2 = 1+floor(Ne*rand);
   v1 = i_srand(e1);
   v2 = j_srand(e1);
   v3 = i_srand(e2);
   v4 = j_srand(e2);
   if (v1~=v3)&&(v1~=v4)&&(v2~=v4)&&(v2~=v3)
      if rand>0.5
        if (srand(v1,v3)==0)&&(srand(v2,v4)==0)
            srand(v1,v2)=0;
            srand(v3,v4)=0;
            srand(v2,v1)=0;
            srand(v4,v3)=0;            
            srand(v1,v3)=1;
            srand(v2,v4)=1;
            srand(v3,v1)=1;
            srand(v4,v2)=1;           
            nrew=nrew+1;          
            i_srand(e1)=v1;
            j_srand(e1)=v3;
            i_srand(e2)=v2;
            j_srand(e2)=v4;
        end
      else
         v5=v3;
         v3=v4;
         v4=v5;
         clear v5;       
         if (srand(v1,v3)==0)&&(srand(v2,v4)==0)
            srand(v1,v2)=0;
            srand(v4,v3)=0;
            srand(v2,v1)=0;
            srand(v3,v4)=0;           
            srand(v1,v3)=1;
            srand(v2,v4)=1;
            srand(v3,v1)=1;
            srand(v4,v2)=1;           
            nrew=nrew+1;           
            i_srand(e1)=v1;
            j_srand(e1)=v3;
            i_srand(e2)=v2;
            j_srand(e2)=v4;
         end      
      end
   end
end