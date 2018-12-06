OBJS=cpg.o main.o options.o

all: cpg

cpg: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(LFLAGS) $(LIBS) $(OBJS)

clean:
	rm -f *.o cpg core *~

.c.o:
	$(CC) -Wall $(CFLAGS) -c $<
