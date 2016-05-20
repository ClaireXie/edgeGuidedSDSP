function rmse=calc_rmse_roi(input, gt)
mask = gt > 0;
sum_value = (input(:)-gt(:)).^2.*mask(:);
rmse = sum(sum_value(:))/sum(mask(:));
rmse = sqrt(rmse);

%rmse=sqrt(mean((input(:)-gt(:)).^2));