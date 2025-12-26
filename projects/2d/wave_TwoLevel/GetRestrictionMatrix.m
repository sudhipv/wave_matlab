% Get Restriction Matrix

function [R] = GetRestrictionMatrix(numpartition,node_gamma,node_gamma_part,level)


for n=1:numpartition
    
     n_g = nonzeros(unique(node_gamma(:,1)));
     
    if(level == 1)
     n_g = nonzeros(unique(node_gamma(:,:,n)));
    end
    
    n_g_p = unique(node_gamma_part(any(node_gamma_part(:,1,n),2),:,n));
    
      
    for i=1:size(n_g,1)
        
        r_ngp = find(n_g(i) == n_g_p);
        
        for j=1:size(n_g_p,1)
            
                if(any(r_ngp) && (j== r_ngp))
                 R(j,i,n) = 1;
                 r_ngp = 0;
                else
                 R(j,i,n) = 0; 
                end
                
        end
        
    end
    
end

end