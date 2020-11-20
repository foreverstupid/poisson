CC=mpixlc			# mpixlc_r   for using OMP on BG/P
CFLAGS=-O2			# -qsmp=omp  for using OMP on BG/P
                    # -qarch=pwr8 for Polus
					# -qarch=450d for BG
LIBS=-lm
LDFLAGS=

SRC=matrix.c solve.c operator.c io.c process.c problem.c
NOSRC=definitions.h
OBJ=$(SRC:%.c=%.o)
NAME=poisson

%.o: %.c %.h $(NOSRC)
	$(CC) $(CFLAGS) -c $< -o $@

$(NAME): main.c $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f $(NAME)