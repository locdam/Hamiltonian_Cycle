function [eidx] = NextEdge(G, S)
    
    [mi, pos] = min(G.Edges.Weight(S));
    eidx = S(pos);
       
end