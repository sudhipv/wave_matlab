%Extract gmsh with partition
% Only for gmsh generated partitoned files - Seperate file for each domain

function [p,e,t, GLnode, npart,point_exact,triangle_exact,node_gamma,...
    node_interior,node_interfacewhole,node_remainingwhole,R,D,...
    node_corner,node_remaining,node_c,node_r,Rc,Rr,Bc,node_gamma_OG,node_interior_OG,node_corner_OG,node_remaining_OG] = ExtractGmshPart_puff_TwoLevel(filename)

file = fullfile(filename);

% 4th line gives number of nodes in the mesh

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
  
  triangle_exact(i,1:3) = triangles(i,partnumtri+2:partnumtri+4);

end
triangle(:,4,j) = partID(j);
end
triangle_exact(:,4) = 1;

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
      p_part(1,n+k:n+2+k,j) = triangles(i,partnumtri+2:partnumtri+4);
      k = k+2;
  end
  
end

p_part(2,:,j) = partID(j);

end

% p_part contains - Node numbers by partition

for j=1:numpartition
    
  p_u = nonzeros(unique(p_part(1,:,j)))';
  
  for i=1:size(p_u,2)   
  point(i,:,j) =  points(p_u(1,i),:) ;
  end
    
end

% point contains node numbers and their co-ordinates in each partition. 

% Seperating out edges and triangles to suit puffin format

p = point;

e = edge;

t = triangle;


% Extracting the interface and remaining nodes and constructing R matrix for each partition 

[node_gamma,node_interior,node_interfacewhole,node_remainingwhole,R,D,node_corner,node_remaining,...
 node_c,node_r,Rc,Rr,Bc,node_gamma_OG,node_interior_OG,node_corner_OG,node_remaining_OG] = GetInterfaceNodes_puff_Twolevel(numpartition,partID,p_part,triangles,edges);


[e_R,t_R,globalnodenum] = GetRenamedEdgeTriangles(numpartition,point,e,t);


e = e_R;
t = t_R;
GLnode = globalnodenum;
npart = numpartition;
point_exact = points(:,2:3)';
triangle_exact = triangle_exact';


end








