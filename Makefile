LIBDIR = /localhome/localuser/lib
INCDIR = /localhome/localuser/include

CC = cc
OFILES = fixalign.o
LIBS   = -lbiop -lgen -lm
CFLAGS = -g -ansi -Wall

fixalign : $(OFILES)
	$(CC) $(CFLAGS) -o $@ $(OFILES) -L $(LIBDIR) $(LIBS)

.c.o :
	$(CC) $(CFLAGS) -c $< -I $(INCDIR)

clean :
	rm $(OFILES)
