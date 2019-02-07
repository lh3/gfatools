CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2
CPPFLAGS=
INCLUDES=	-I.
OBJS=		gfa-io.o gfa-asm.o gfa-util.o
PROG=		gfaview
LIBS=		-lz

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

gfaview:$(OBJS) gfaview.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

gfa-chk:gfa-chk.l
		lex $< && $(CC) -O2 lex.yy.c -o $@

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session* gfa-chk

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE
