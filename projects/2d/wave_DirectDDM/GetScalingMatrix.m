function [D] = GetScalingMatrix(numpartition,node_interface,node_interface_part)

nodes = unique(nonzeros(node_interface(:,1)));

for i=1:length(nodes)

           nodes(i,2) = length(unique(node_interface(find(node_interface(:,1) == nodes(i)),2:3)));

end

for n=1:numpartition
    
nodes_part = unique(nonzeros(node_interface_part(:,1,n)));

for i = 1:length(nodes_part)
   
    D(i,i,n) = 1/nodes(find(nodes(:,1)== nodes_part(i)),2); 
    
    
end 

end

end