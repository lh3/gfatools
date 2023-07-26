CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -std=c99 -O2
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kalloc.o gfa-base.o gfa-io.o gfa-aug.o gfa-sub.o gfa-asm.o gfa-util.o \
			gfa-bbl.o gfa-ed.o gfa-sql.o
EXE=		gfatools
LIBS=		-lz

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -lpthread -ldl -lm
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

gfatools:main.o sys.o libgfa1.a
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

libgfa1.a:$(OBJS)
		$(AR) -csr $@ $(OBJS)

gfa-chk:gfa-chk.l
		lex $< && $(CC) -O2 lex.yy.c -o $@

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM session* gfa-chk

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

gfa-asm.o: gfa-priv.h gfa.h kvec.h kdq.h kalloc.h
gfa-aug.o: gfa-priv.h gfa.h ksort.h
gfa-base.o: gfa-priv.h gfa.h kstring.h khashl.h kalloc.h ksort.h
gfa-bbl.o: gfa-priv.h gfa.h kalloc.h ksort.h kvec.h
gfa-ed.o: gfa-priv.h gfa.h kalloc.h ksort.h khashl.h kdq.h kvec-km.h
gfa-io.o: kstring.h gfa-priv.h gfa.h kseq.h
gfa-sql.o: kstring.h gfa-priv.h gfa.h
gfa-sub.o: gfa-priv.h gfa.h kalloc.h kavl.h khashl.h ksort.h kvec.h
gfa-util.o: gfa-priv.h gfa.h kvec.h ksort.h kdq.h kalloc.h khashl.h
kalloc.o: kalloc.h
main.o: ketopt.h gfa-priv.h gfa.h kalloc.h kseq.h
