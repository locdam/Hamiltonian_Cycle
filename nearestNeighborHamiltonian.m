function [L,wt] = nearestNeighborHamiltonian(G)

% check if vertices have names
if (~sum(ismember(G.Nodes.Properties.VariableNames,'Name')))
    % if not, give names using its indices
    Vnames = int2str(1:numnodes(G));
    G.Nodes.Name = split(Vnames);
end

% check if edges have names
if (~sum(ismember(G.Edges.Properties.VariableNames,'Name')))
    % if not, give names using its indices
    Enames = int2str(1:numedges(G));
    G.Edges.Name = split(Enames);
end

v_id = 1;
% set the dfnumber of all vertices to -inf
G.Nodes.dfN = -inf(numnodes(G),1);

% Let T = the vertex with id v_id
T = graph;
T = addnode(T,1);

% record the original id for the vertex in G in the origId attribute of the
% nodes of T
T.Nodes.origId(1) = v_id;
T.Nodes.Name(1) = G.Nodes.Name(v_id);


% initiate the counting of dfnumber
currentDf = 0;

% set the dfnumber for the starting vertex in G and in T
G.Nodes.dfN(v_id) = currentDf;
T.Nodes.dfN(1) = currentDf;

% the first set of frontier edges are the edges from the vertex with id
% v_id
[S,nV] = outedges(G,v_id);

S = S(nV~=v_id);

while ~isempty(S)
    currentDf = currentDf+1;
    
    % edge idx (in G) of the next edge
    eidx = NextEdge(G,S);

%     endpints of the chosen next edge
    endpts = G.Edges.EndNodes(eidx,:);
    endpts = findnode(G,{endpts{1} endpts{2}});
    
    % next vertex and its tree node
    if (~isinf(G.Nodes.dfN(endpts(1)))) % endpts(1) is a tree node
        w_id = endpts(2);
        pre_id = endpts(1);
    else % endpts(2) is a tree node
        w_id = endpts(1);
        pre_id = endpts(2);
    end
    % add the new node to the tree
   
     % add the new node to the tree
    newNode = table(G.Nodes.Name(w_id), w_id, currentDf,'VariableNames', {'Name','origId', 'dfN'});
    T = addnode(T,newNode);
    
    % create the edge and its attributes (endpts and original id in G) to be added in T
    newEdge = table([G.Nodes.dfN(pre_id)+1,currentDf+1],G.Edges.Name(eidx),eidx,G.Edges.Weight(eidx),'VariableNames', {'EndNodes','Name','origId','Weight'});
    T = addedge(T,newEdge);
    % record the dfnumber in G and T for the next discovered vertex
    G.Nodes.dfN(w_id) = currentDf;
    
    % update the set of frontier edges
    S = updateFrontierEdge(G,S);
    
end
    eidx = findedge(G, T.Nodes.origId(1), T.Nodes.origId(end));
   
    
    % create the edge and its attributes (endpts and original id in G) to be added in T
    newEdge = table([G.Nodes.dfN(pre_id)+1,currentDf+1],G.Edges.Name(eidx),eidx,G.Edges.Weight(eidx),'VariableNames', {'EndNodes','Name','origId','Weight'});
    T = addedge(T,newEdge);
    % record the dfnumber in G and T for the next discovered vertex
    G.Nodes.dfN(w_id) = currentDf;
    [~,asort]=sort(T.Nodes.origId); 
    b= T.Nodes.dfN;
    L=b(asort);

    wt = sum(T.Edges.Weight);
end