function B = nmodeproduct(A,M,n)
dimvec = size(A);
n = fix(n);
if (length(dimvec)<n || n<1)
    error('nmodeproduct: n is not within the order range of tensor A ');
end
if (size(M,2) ~= dimvec(n))
    error('nmodeproduct: dimension n of tensor A is not equal to dimension 2 of matrix M');
end
% shift A to prepare flattening: (i.e. make dimension 1 (columns) to 'n', the one we would like to replace)
Ash = shiftdim(A,n-1);
% save the target dimensions of B (we replace the 1st dimension, because
% thats the one affected by the matrix multiplication 
% i.e. this dimension changes from I_n to J
dimvecB = size(Ash);
dimvecB(1) = size(M,1);
% multiply while flattening.. i.e. we first flatten the matrix, so that we
% have a matrix; as an array of flattened vectors drawn from the tensor
% second we multiply those vectors with our matrix, resulting in a
% dimension change of the output vector. the output vectors are then again
% saved as a matrix, representing an array of those vectors.
B = M*Ash(:,:);
% wrap the flattened vector-array back into the previously saved tensor shape
B = reshape(B,dimvecB);
% shift the dimensions back! so that only dimension n has changed from I_n to J
B = shiftdim(B,length(dimvecB)-n+1);
% and we re done!