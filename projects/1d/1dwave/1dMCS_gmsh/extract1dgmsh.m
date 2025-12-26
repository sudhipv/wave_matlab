

%%% extract 1D mesh


% ExtractGmsh
file = fullfile('./axialbar.msh');

% 4th line gives number of nodes in the mesh

nodenumber = dlmread(file,' ',[4 0, 4 0]);

% Lines from 0-4 in .msh file gives descriptions

nodes = dlmread(file,' ',[5 0 nodenumber+4 3]);

% 3 Lines after node number are descriptions
% 3rd line after nodes provides number of elements

elenumber = dlmread(file,' ',[nodenumber+7 0, nodenumber+7 0]);

elements = dlmread(file,' ',[nodenumber+8 0 elenumber+nodenumber+7 6]);

% Here line elements will have last digit assigned 0 by code which is not
% required

% extracting points/nodes edges and traingles from the extracted arrays
% from gmsh

% For 2D only 2 and 3rd columns which gives x and y co-ordinates are taken

points = nodes(1:nodenumber,2:3); % <x-co-ordinate, y co-ordinate>

elements = elements(1:elenumber,6:7);
p = points;
e = elements;

%Save Mesh data
dlmwrite('points.txt', p,' ')
dlmwrite('elements.txt', e,' ')
dlmwrite('nodes.txt',nodes(:,1),' ')

