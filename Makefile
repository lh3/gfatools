CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kalloc.o gfa-base.o gfa-io.o gfa-aug.o gfa-sub.o gfa-asm.o gfa-util.o
PROG=		gfatools
LIBS=		-lz

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

gfatools:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

gfa-chk:gfa-chk.l
		lex $< && $(CC) -O2 lex.yy.c -o $@

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session* gfa-chk

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE
