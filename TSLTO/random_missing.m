function [sparse_tensor,Omega] = random_missing(T,missing_rate)
sizeT = size(T);
random_tensor = rand(sizeT);
R = round(random_tensor + 0.5 - missing_rate);
sparse_tensor = T .* R;
Omega = (R > 0);