CC=gcc
CFLAGS= -O2

progs = box

headers = 

obj = viewfactors.o

exes = $(addsuffix .exe, $(progs))

all: $(progs) $(obj)
$(obj): $(headers)
$(progs): $(obj) $(headers)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

.c:
	$(CC) $(CFLAGS) -o $@ $*.c $(obj) -lm

.c.exe:
	$(CC) $(CFLAGS) -o $@ $*.c $(obj) -lm

.f.o:
	$(FPP) $< $*F.f
	$(F77) -c $(F77FLAGS) -o $@ $*F.f
	rm $*F.f
.f:
	$(FPP) $(DEFFLAGS) -I$(incdir) $< $*F.f
	$(F77) $(F77FLAGS) -o $@ $*F.f $(libjsl) $(teclib)
	rm $*F.f
	cp $@$(EXT) $(bindir)

clean:
	rm -f $(progs) $(obj)
