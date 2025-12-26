% Function to rename edges and triangles for assembling matrix domain level

function [e,t,globalnodenum] = GetRenamedEdgeTriangles(numpartition,point,e,t)

num = 0;

 for num = 1:numpartition
     
     pt = nonzeros(point(:,1,num));
     
     for n = 1:length(pt)

      [r_e,c_e] = ind2sub(size(e(:,:,num)),find(e(:,:,num) == pt(n)));
      [r_t,c_t] = ind2sub(size(t(:,:,num)),find(t(:,:,num) == pt(n)));

              for i=1:length(r_e)

                  if(c_e(i) ~= 1 &&c_e(i) ~=4)
                    e(r_e(i),c_e(i),num) =  n; 
                  end

              end

              for j=1:length(r_t)

                    if(c_t(j) ~= 4)
                     t(r_t(j),c_t(j),num) =  n; 
                    end

              end

     end
     
 pt = point(:,1,num);  % Assigning pt back to get normal size containing zeros 
 globalnodenum(:,:,num) = pt;
 
 end
 
end