function PrintFilter(file, f, in, out, h)
    fprintf(file,'FILTER Filter%03d Length=%d In=%d Out=%d \n ',f,size(h,1),in,out);
    fprintf(file,[repmat('%.5f ',1, 20) '\n '],h(1:end));
    