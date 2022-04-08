function [RMSE,maxerr] = rmse(err)
    [m,n] = size(err);
    if m>n  
        RMSE = sqrt(1/length(err)*sum(err.^2,1));  %返回行向量
    else
        RMSE = sqrt(1/length(err)*sum(err.^2,2)); %返回列向量
    end
    
    maxerr = max(abs(err));
end