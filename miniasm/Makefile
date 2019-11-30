CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -std=c99 -O2
CPPFLAGS=
INCLUDES=	-I..
OBJS=		misc.o paf.o hit.o utg.o
EXE=		miniasm
LIBS=		-lz

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

miniasm:$(OBJS) ../sys.o main.o ../libgfa1.a
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

../libgfa1.a:
		(cd .. && $(MAKE))

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE
