function rmse=calc_rmse(input, gt)
rmse=sqrt(mean((input(:)-gt(:)).^2));