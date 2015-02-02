function errPer=calc_err(input,gt, th)
tmp=input-gt;
errPer=100*sum(sum((abs(tmp)>th)))/(size(tmp,1)*size(tmp,2));