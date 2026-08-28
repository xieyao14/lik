% continuous influence kernel
function eff = truf_test141(k1, k2, tau, lag)
    %%
    K = 5;
    edges = [2,1;1,2;2,3;5,4;1,5;4,3;2,2;5,5;];
    adj = zeros(K,K);
    for i = 1:size(edges,1)
        adj(edges(i,1),edges(i,2))=1;
    end
    G = digraph(adj);

    peak = load("peak.mat");
    freq = load("freq.mat");
    peak = peak.peak;
    freq = freq.freq;

    eff = 0;
    for i = 1:size(G.Edges.EndNodes,1)
        if k2==G.Edges.EndNodes(i,1) && k1==G.Edges.EndNodes(i,2)  
            eff = (cos(freq(i).*(tau+2))+0.75)*0.35.*exp(-20.*(lag-peak(i)).^2);
        end
    end

end