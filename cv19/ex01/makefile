CC = gcc

EXE = hough_transform gradient_map
OBJ = hough_transform.o gradient_map.o
SRC = gradient_map.c hough_transform.c 

OPTIMIZATION = -O3

all: $(EXE) #gradient_map hough_transform

grad:
	./gradient_map < ./in_grad.txt
	eog coins-grad.pgm

hough:
	./hough_transform < ./in_hough.txt
	eog coins-circles.pgm

%.o: %.c
	$(CC) $< -o $@ -c $(OPTIMIZATION) $(DBGCFLAGS) -lm

gradient_map: gradient_map.o
	$(CC) $^ -o $@ $(OPTIMIZATION) $(DBGCFLAGS) -lm

hough_transform: hough_transform.o
	$(CC) $^ -o $@ $(OPTIMIZATION) -lm

debug: DBGCFLAGS = -DDEBUG -g
debug: OPTIMIZATION = -O0			# optimization level = 0 (otherwise breakpoints will not work with optimized code)
debug: $(OBJ) $(EXE)

gif:
	convert -delay 100 -loop 0 coins.pgm coins-grad.pgm coins-circles.pgm coins-circles-only.pgm  coins.gif

clean:	
	rm -f *.o
	rm -f *.a
	rm -f $(OBJ) $(EXE)

1c:
	rm -f $(OBJ) $(EXE)

.PHONY: all clean lc debug