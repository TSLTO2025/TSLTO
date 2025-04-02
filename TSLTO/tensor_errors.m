function [rmse, mape, mae] = tensor_errors(A, B, Omega)
    % 输入:
    % A, B: 两个大小相同的张量
    % 输出:
    % rmse: 均方根误差 (Root Mean Square Error)
    % mape: 平均绝对百分比误差 (Mean Absolute Percentage Error)
    
    % 检查A和B尺寸是否相同
    if ~isequal(size(A), size(B))
        error('A and B must have the same dimensions.');
    end
    
    %RMSE
    diff = A - B;%算完整张量
    %diff = A(Omega==0) - B(Omega==0);%只算缺失部分
    rmse = sqrt(mean(diff(:).^2));
    
    %MAPE
    if any(A(:) == 0)
        warning('A contains zero elements, MAPE may be undefined for these values.');
        %避免除0
        non_zero_indices = A ~= 0;
        mape = mean(abs((A(non_zero_indices) - B(non_zero_indices)) ./ A(non_zero_indices))) * 100;
    else
        mape = mean(abs((A(:) - B(:)) ./ A(:))) * 100;
    end

    %MAE
    absolute_error = abs(A - B);
    mae = mean(absolute_error(:)); 
    
    fprintf('RMSE: %.4f\n', rmse);
    fprintf('MAPE: %.4f%%\n', mape);
    fprintf('MAE: %.4f\n', mae);
end
