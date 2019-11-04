CC = gcc

EXE = hough_transform gradient_map
OBJ = hough_transform.o gradient_map.o
SRC = gradient_map.c hough_transform.c 

all: $(EXE) #gradient_map hough_transform

grad:
	./gradient_map < ./in_grad.txt
	display coins-grad.pgm

hough:
	./hough_transform < ./in_hough.txt
	display coins-class-map.pgm

%.o: %.c
	$(CC) $< -o $@ -c -O3 -lm

gradient_map: gradient_map.o
	$(CC) $^ -o $@ -O3 -lm

hough_transform: hough_transform.o
	$(CC) $^ -o $@ -O3 -lm

clean:	
	rm -f *.o
	rm -f *.a
	rm -f $(OBJ) $(EXE)

1c:
	rm -f $(OBJ) $(EXE)

.PHONY: all clean lc debug