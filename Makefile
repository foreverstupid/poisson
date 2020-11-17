CC=mpicc
CFLAGS=-Wall -O5
LIBS=-lm
LDFLAGS=

SRC=matrix.c solve.c operator.c output.c process.c problem.c
NOSRC=definitions.h
OBJ=$(SRC:%.c=%.o)
NAME=poisson

%.o: %.c %.h $(NOSRC)
	$(CC) $(CFLAGS) -c $< -o $@

$(NAME): main.c $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRC)
	$(CC) -MM $^ > $@

clean:
	rm -f deps.mk
	rm -f *.o
	rm -f $(NAME)
