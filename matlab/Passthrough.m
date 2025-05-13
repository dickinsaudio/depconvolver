
filter=0;
file = fopen('PASSTHROUGH.bin','w');


for (i=1:64)
            PrintFilter2(file,filter,i,i,1,true);
            filter = filter+1;
end;
fclose(file);
