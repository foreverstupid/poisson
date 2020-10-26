CC=mpicc
CFLAGS=-Wall -g -fno-exceptions
LIBS=-lm
LDFLAGS=

SRC=matrix.c solve.c operator.c
NOSRC=definitions.h
OBJ=$(SRC:%.c=%.o)
NAME=poisson

%.o: %.c %.h $(NOSRC)
	$(CC) $(CFLAGS) -c $< -o $@

$(NAME): main.c $(OBJ) $(NOSRC)
	$(CC) $(CFLAGS) $(INCDIR) $(LDFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f $(NAME)
