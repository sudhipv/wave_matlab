 
% Finding corner and remianing nodes from interface nodes


function [node_corner,node_remaining,node_corner_part,node_remaining_part] = GetCoarseGrid(numpartition,node_interface,node_interface_part,node_cor,edges)

 node_inter = nonzeros(unique(node_interface(:,1)));
 node_corner = nonzeros(unique(node_cor(:,1)));
 n = 1;
 
 for i = 1:length(node_inter)
     
 if(any(node_inter(i)==edges(:,4:5)))
    node_corner(length(node_corner)+1) = node_inter(i);
 else
     if(~(any(node_inter(i)== node_corner(:,1))))
     node_remaining(n) = node_inter(i);
     n = n+1;
     end
 end
 
 end
 
  node_remaining = node_remaining';
  
  
  node_corner = nonzeros(unique(node_corner));
  
 % Making Partitioned Corner and Remaining Nodes
  
 for j = 1: numpartition
    
   node_inter_part = nonzeros(unique(node_interface_part(:,j)));
   nc = 1;
   nr = 1;
   
           for i = 1:length(node_inter_part)

                 if(any(node_inter_part(i)==node_corner(:,1)))
                    node_corner_part(nc,:,j) = node_inter_part(i);
                    nc = nc+1;
                 else
                     node_remaining_part(nr,:,j) = node_inter_part(i);
                     nr = nr+1;
                 end
          end
 
 end
 
  
end