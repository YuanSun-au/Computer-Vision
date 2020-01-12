CC = gcc
EXE = ced
OBJ = ced64.o
SRC = diff_tensor.c

OPTIMIZATION = -O3

$(EXE): $(SRC)
	$(CC) $(OPTIMIZATION) -o $@ $(OBJ) $< -lm -fPIC
	# gcc -O3 -o ced ced64.o diff_tensor.c -lm

debug: DBGCFLAGS = -DDEBUG -g
debug: OPTIMIZATION = -O0			# optimization level = 0 (otherwise breakpoints will not work with optimized code)
debug: $(EXE)


clean:
	rm -f $(EXE)