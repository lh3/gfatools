CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -std=c99 -O2
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kalloc.o gfa-base.o gfa-io.o gfa-aug.o gfa-sub.o gfa-asm.o gfa-util.o
PROG=		gfatools
LIBS=		-lz

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

gfatools:$(OBJS) main.o sys.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

gfa-chk:gfa-chk.l
		lex $< && $(CC) -O2 lex.yy.c -o $@

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session* gfa-chk

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

gfa-asm.o: gfa-priv.h gfa.h kvec.h kdq.h
gfa-aug.o: gfa-priv.h gfa.h ksort.h
gfa-base.o: gfa-priv.h gfa.h khash.h kalloc.h ksort.h
gfa-io.o: kstring.h gfa-priv.h gfa.h kseq.h
gfa-sub.o: gfa.h kalloc.h kavl.h khash.h ksort.h
gfa-util.o: gfa.h kvec.h
kalloc.o: kalloc.h
main.o: ketopt.h gfa.h kseq.h
