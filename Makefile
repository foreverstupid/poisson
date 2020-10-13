CC=gcc
CFLAGS=-Wall -g -fno-exceptions
LIBS=-lm
LDFLAGS=

SRC=vector.c
OBJ=$(SRC:%.c=%.o)
NAME=poisson

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

$(NAME): main.c $(OBJ)
	$(CC) $(CFLAGS) $(INCDIR) $(LDFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f $(NAME)
