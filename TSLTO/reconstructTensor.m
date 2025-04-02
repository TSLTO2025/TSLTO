function x = reconstructTensor(G,U)
    x = double(G);
    for i = 1:3
        x = nmodeproduct(x,U{i},i);
    end
end