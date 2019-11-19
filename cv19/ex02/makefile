CC = gcc

EXE = lkTemplate
OBJ = lkTemplate.o
SRC = lkTemplate.c 

OPTIMIZATION = -O3

all: $(EXE)

pig:
	./lkTemplate  < ./pig-params.txt
	eog pig-flow-class.pgm
	eog pig-flow-mag.pgm

sphere:
	./lkTemplate  < ./sphere-params.txt
	eog sphere-flow-class.pgm
	eog sphere-flow-mag.pgm

%.o: %.c
	$(CC) $< -o $@ -c $(OPTIMIZATION) $(DBGCFLAGS) -lm

$(EXE): $(OBJ)
	$(CC) $^ -o $@ $(OPTIMIZATION) $(DBGCFLAGS) -lm


debug: DBGCFLAGS = -DDEBUG -g
debug: OPTIMIZATION = -O0			# optimization level = 0 (otherwise breakpoints will not work with optimized code)
debug: $(OBJ) $(EXE)

# gif:
# 	convert -delay 100 -loop 0 coins.pgm coins-grad.pgm coins-circles.pgm coins-circles-only.pgm  coins.gif

clean:	
	rm -f *.o
	rm -f *.a
	rm -f $(OBJ) $(EXE)

1c:
	rm -f $(OBJ) $(EXE)

.PHONY: all clean lc debug