
filter=0;
file = fopen(sprintf('PASSTHROUGH.txt',set),'wt');


for (i=1:64)
            PrintFilter2(file,filter,i,i,1);
            fprintf(file,'\n\n');
            filter = filter+1;
end;
fclose(file);
