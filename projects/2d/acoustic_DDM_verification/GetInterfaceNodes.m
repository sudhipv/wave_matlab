  

%Extract gmsh with partition
clear all
file = fullfile('squarepart4.msh');

nodenumber = dlmread(file,' ',[4 0, 4 0]);
nodes = dlmread(file,' ',[5 0 nodenumber+4 3]);
% 3 Lines after node number are descriptions
% 3rd line after nodes provides number of elements

elenumber = dlmread(file,' ',[nodenumber+7 0, nodenumber+7 0]);
maxtagnum = max(dlmread(file,' ',[nodenumber+8 2 elenumber+nodenumber+7 2]));
elements = dlmread(file,' ',[nodenumber+8 0 elenumber+nodenumber+7 maxtagnum+5]);

% Here line elements will have last digit assigned 0 by code which is not
% required

points = nodes(1:nodenumber,1:3); % <node number, x-co-ordinate, y co-ordinate>

% Loop to find the position of starting of triangles from element array
for i=1:elenumber
    
    if(elements(i,2) == 2)
        tristart = i;
        break
    end
end


% Handling Partitioned MSH file

% tag number for all edges will be same, they belong to only 

edges = elements(1:tristart-1,5:9); % <geometrical tag ,number of part, partID, node 1, node 2 > 

triangles = elements(tristart:elenumber, 6:maxtagnum+6); % <partnumber, partID, node1, node 2, node3>

% 6th tag is number of partitions the element is shared with

% 3rd tag is the number of tags, First two tags are physical and
% geometrical identities

partID = unique(elements(tristart:elenumber,7));
numpartition = length(partID);


% Finding out the size of triangle and assigning it to zero to avaoid using
% Dynamic Memory Allocation

% Triangles, edges, points are the matrices with all the triangles,edges and points in domain
%where as triangle,edge and point are 3D matrices with seperated out
%quantities in terms of patitions

for j = 1:numpartition
    n = 0;
    ns = 0;
   
for i=1:size(triangles,1)
    
  if(triangles(i,2) == partID(j))
      n = n +1;
  end
  
  
end
num(j) = n;
triangle = zeros(num(j),4,j);

end

% Seperated out the triangles with respect to partition number in a 3D
% matrix with 3rd indice corresponding to partition it belongs traingle(:,:,1) - 1st partition 

% traingle(:,:,1) < node1, node2, node3> in partition 1
% traingle(:,:,2) < node1, node2, node3> in partition 2 and so on

for j=1:numpartition
    n = 0;
    ns = 0;
    ni = 0;
for i=1:size(triangles,1)
    
    partnumtri = triangles(i,1);
    
  if(triangles(i,2) == partID(j))
      n = n +1;
      triangle(n,1:3,j) = triangles(i,partnumtri+2:partnumtri+4);
  end
  

end
triangle(:,4,j) = partID(j);
end

% Seperating out the edges

for j= 1:numpartition
    n = 0;
    
for i = 1:size(edges,1)
    
    if(edges(i,3) == partID(j))
      n = n +1;
      edge(n,1,j) = edges(i,1);
      edge(n,2:3,j) = edges(i,4:5);
    end
    
end
    edge(:,4,j) = partID(j);
end

% edge(:,:,1) = <physical edge number, node 1, node 2> partition 1
% edge(:,:,2) = <physical edge number, node 1, node 2> partition 2


% Seperating out the points

for j=1:numpartition
    n = 0;
    k = 0;
for i=1:size(triangles,1)
    
    partnumtri = triangles(i,1);
    
  if(triangles(i,2) == partID(j))
      n = n +1;
      p(1,n+k:n+2+k,j) = triangles(i,partnumtri+2:partnumtri+4);
      k = k+2;
  end
end
p(2,:,j) = partID(j);
end


% point contains just node numbers in each partition. Since the shared
% nodes are not seperated out it contains repetitive node numbers.
msizepoint = 0;
for j=1:numpartition
    
  p_u = nonzeros(unique(p(1,:,j)))';
  
  if(size(p_u,2) > msizepoint)
  msizepoint = size(p_u,2);
  end

  point(:,:,j) = zeros(msizepoint,3,1);
  
  for i=1:size(p_u,2)   
  point(i,:,j) =  points(p_u(1,i),:) ;
  end
    
end

% Extracting Interface nodes

% Finding out shared triangles
ns = 0;
for i=1:size(triangles,1)
    
      
  if(triangles(i,1) ~= 1)
      ns = ns +1;
      triangle_shared(ns,:) = triangles(i,:);
  end
  
    
end

% Extracting shared triangles by partition ID
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
% 2. Check for the partition ID the triangle is shared with
% 3. Now search for the same node in that part ID in whole shared triangles
% 4. If its possible to find out the same node atleast one time, its an
% interface node

% Corner Nodes are also found from the same loop
% If the partition number is more than 2 one of the node is a corner node.
% If the node is an interface node, check for the 2 nd partition ID and
% match the corresponding node in that part ID to confirm whether 
%its a corner node or not 

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


% node_interface = unique(node_interface);
% node_corner = unique(node_corner);

dlmwrite('node_interface.txt', node_interface,' ')
dlmwrite('node_corner.txt', node_corner,' ')
%dlmwrite('triangleshared_part.txt', triangleshared_part,' ')

for n=1:numpartition
  
    node_part1 = find(node_interface(:,2) ==partID(n));
    node_part2 = find(node_interface(:,3) ==partID(n));
    node_part = unique(vertcat(node_part1,node_part2));
    for i=1:size(node_part,1)
    node_interface_part(i,:,n)=node_interface(node_part(i),1);
    end
end

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






