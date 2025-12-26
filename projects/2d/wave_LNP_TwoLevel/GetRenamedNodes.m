
% Function to rename nodes for assembling matrix domain level

function [nodes_1R,nodes_2R] = GetRenamedNodes(numpartition,p_part,nodes_1,nodes_2)

 for num = 1:numpartition
     
  p_u = nonzeros(unique(p_part(1,:,num)))';
     
 for n = 1:size(p_u,2)
     
  [r_r,c_r] = ind2sub(size(nodes_1(:,:,num)),find(nodes_1(:,:,num) == p_u(n)));
  [r_i,c_i] = ind2sub(size(nodes_2(:,1,num)),find(nodes_2(:,1,num) == p_u(n)));
    
          for i=1:size(r_r,1)

              nodes_1R(r_r(i),c_r(i),num) =  n; 

          end
  
          for i=1:size(r_i,1)

              nodes_2R(r_i(i),c_i(i),num) =  n; 

          end
   
 end
 
 end

end