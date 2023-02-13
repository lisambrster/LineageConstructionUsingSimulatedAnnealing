function ggboth = CombineGraphs(gg1,gg2)

% combine two graphs (gg1,gg2) into one graph

% figure;
% plot(gg1,'layout','layered');  
% figure;
% plot(gg2,'layout','layered');  


ggboth = graph;
nNodes = size(gg1.Nodes);
if (nNodes == 0)
    ggboth = gg2;
else
    for iNode = 1:nNodes
        ggboth = ggboth.addnode(gg1.Nodes{iNode,1});
    end
    
    nNodes = size(gg2.Nodes);
    for iNode = 1:nNodes
        node = gg2.Nodes{iNode,1};
        % check that not already in new graph
        k = findnode(ggboth,node);
        if (k < 1)
            ggboth = ggboth.addnode(node);
        end
    end
    
    
    nEdges = size(gg1.Edges);
    for iEdge = 1:nEdges
        node1 = gg1.Edges{iEdge,1}(1,1);
        node2 = gg1.Edges{iEdge,1}(1,2);
        ggboth = ggboth.addedge(node1,node2);
    end
    
    nEdges = size(gg2.Edges);
    for iEdge = 1:nEdges
        node1 = gg2.Edges{iEdge,1}(1,1);
        node2 = gg2.Edges{iEdge,1}(1,2);
        ggboth = ggboth.addedge(node1,node2);
    end
    
    % show each of the new graphs
    %figure;
    %plot(ggboth,'layout','layered');  
   
  end
end

