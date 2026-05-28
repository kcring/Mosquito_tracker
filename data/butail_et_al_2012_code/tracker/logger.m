function logger(logfile, str)

lf=fopen(logfile, 'a');

fprintf(lf, '\n[%s] %s', datestr(now, 'yy-mmm-dd HH.MM.SS'), str);

fclose(lf);

