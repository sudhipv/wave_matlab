  
% Extracting Interface nodes

function [node_gamma,node_interior,node_interfacewhole,node_remainingwhole,R]... 
= GetInterfaceNodes_puff(numpartition,partID,p_part,triangles)

% Finding out shared triangles

ns = 0;

for i=1:size(triangles,1)
    
      
  if(triangles(i,1) ~= 1)
      ns = ns +1;
      triangle_shared(ns,:) = triangles(i,:);
  end
  
    
end

% Seperating shared triangles by partition ID

for j=1:numpartition
    
    n= 0;
   
    for ns = 1:size(triangle_shared,1)

        partnumtri_share = triangle_shared(ns,1);

              if(triangle_shared(ns,2) == j)

                  n = n +1;
                  triangleshared_part(n,:,j) = triangle_shared(ns,partnumtri_share+2:partnumtri_share+4);

              end
              
    end
    
end

%Finding out the interface nodes by searching through shared triangles.
%The logic is as follows :
% 1. Take any particular node of a shared triangle
% 2. Get the partition ID's of the triangle to which it is shared. One with
% - ve signs
% 3. Now search for the same node in that part IDs in whole shared triangles
% 4. If its possible to find out the same node atleast one time, its an
% interface node
% 5. This method ensures that every interface node is captured even though
% it repeats the same nodes many times.

% Corner Nodes are also found from the same loop
% If the partition ID's are more than 2, one of the node is a corner node.
% if its possible to find the same node in 2 or more than 2 part IDs of the
% shared triangle part ID's then its a corner node.

% Node at the physical edge is not being picked in here and thus not
% included as corner node.

ni = 0;
nc = 0;
for ns = 1:size(triangle_shared,1)
  
    ns_partnum = triangle_shared(ns,1);
    ns_partID = -(triangle_shared(ns,3:ns_partnum+1));
    ns_trueID = triangle_shared(ns,2);
    ns_search(1:3) = triangle_shared(ns,ns_partnum+2:ns_partnum+4);
     interfacenum =0;
     cornernum = 0;
     
   if(ns_partnum >2)
       cornernum = 1;
       ns_partIDCorner = -(triangle_shared(ns,3:ns_partnum+1));
   end
   
   for i = 1:3
       foundInter =0;
       foundCorner = 0;
       
             
 for j=1:size(triangle_shared,1)
       
      nj_partnum = triangle_shared(j,1);
      nj_partID = triangle_shared(j,2);
      
            
      if(any(nj_partID == ns_partID) && (ns_search(i) == triangle_shared(j,nj_partnum+2) ||ns_search(i) == triangle_shared(j,nj_partnum+3) ||...
      ns_search(i) == triangle_shared(j,nj_partnum+4)))
            foundInter = 1;
      end
      
      if(cornernum==1 && any(nj_partID == ns_partIDCorner) && (ns_search(i) == triangle_shared(j,nj_partnum+2) ||ns_search(i) == triangle_shared(j,nj_partnum+3) ||...
      ns_search(i) == triangle_shared(j,nj_partnum+4)))
            foundCorner = foundCorner + 1;
            if(foundCorner ==1)
            nj_partcornerfound = nj_partID;
            end
      end
      
      if(foundInter > 0)
          
          ni = ni+1;
          node_interface(ni,1) = ns_search(i);
          node_interface(ni,2) = nj_partID;
          node_interface(ni,3) = ns_trueID;
          foundInter = 0;
              if(cornernum ==1)
                  if(foundCorner > 1 && (nj_partcornerfound~=nj_partID))
                     nc = nc+1;
                     node_corner(nc,1) = ns_search(i);
                     foundCorner = 0;
                  end
              end
         
      end
        
 end
 end
end


for n=1:numpartition
  
    node_part1 = find(node_interface(:,2) ==partID(n));
    node_part2 = find(node_interface(:,3) ==partID(n));
    node_part = unique(vertcat(node_part1,node_part2));
    for i=1:size(node_part,1)
    node_interface_part(i,:,n)=node_interface(node_part(i),1);
    end
end

p = p_part;

for j=1:numpartition
    r = 0;
  p_u = nonzeros(unique(p(1,:,j)))';
    
  inter_u = nonzeros(unique(node_interface_part(:,1,j)))';
  
  for i=1:length(p_u)
      
      if(~(any(p_u(i)==inter_u)))
          r = r+1;
          node_r(r,:,j) = p_u(i) ;
      end
      
  end
 
end

node_interior = node_r;
node_gamma = node_interface_part;
node_remainingwhole = node_interior;

% Getting Restriction Matrix before renaming

[R] = GetRestrictionMatrix(numpartition,node_interface,node_interface_part);


%renaming nodes

num = 0;

 for num = 1:numpartition
     
  p_u = nonzeros(unique(p(1,:,num)))';
     
 for n = 1:size(p_u,2)
     
  [r_r,c_r] = ind2sub(size(node_r(:,:,num)),find(node_r(:,:,num) == p_u(n)));
  [r_i,c_i] = ind2sub(size(node_interface_part(:,1,num)),find(node_interface_part(:,1,num) == p_u(n)));
    
          for i=1:size(r_r,1)

              node_interior(r_r(i),c_r(i),num) =  n; 

          end
  
          for i=1:size(r_i,1)

              node_gamma(r_i(i),c_i(i),num) =  n; 

          end
   
 end
 
 end

 node_interfacewhole  = nonzeros(unique(node_interface(:,1)));
  
end




