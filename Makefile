CC=mpicc
CFLAGS=-Wall -O3 -Ofast
LIBS=-lm
LDFLAGS=

SRC=matrix.c solve.c operator.c output.c
NOSRC=definitions.h
OBJ=$(SRC:%.c=%.o)
NAME=poisson

%.o: %.c %.h $(NOSRC)
	$(CC) $(CFLAGS) -c $< -o $@

$(NAME): main.c $(OBJ) $(NOSRC)
	$(CC) $(CFLAGS) $(INCDIR) $(LDFLAGS) $^ $(LIBS) -o $@

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRC)
	$(CC) -MM $^ > $@

clean:
	rm deps.mk
	rm -f *.o
	rm -f $(NAME)
