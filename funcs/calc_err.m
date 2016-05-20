function errPer=calc_err(input,gt, th)
msk = (gt>0);
tmp = (input-gt).*msk;
errPer=100*sum(sum((abs(tmp)>th)))/(sum(msk(:)));