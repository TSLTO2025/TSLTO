function mat=ten2mat(ten,dim,k)
mat=reshape(permute(ten,[k,1:k-1,k+1:length(dim)]),dim(k),[]);
end