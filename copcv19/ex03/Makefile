CC = g++

NAME1 = frontend

BIN_NAME1 = $(NAME1)
OBJ_NAME1 = $(NAME1).o
SRC_NAME1 = $(NAME1).cc

NORMAL_CFLAGS = -O3 -funroll-loops -w
#EXTRA_CFLAGS = -fomit-frame-pointer #-static


INCL_PATH =

LIBS_PATH =
LIBS = #-lm

OGL_LIBS = -lGL -lGLU -lglut
OGL_LIBS_PATH = -L/usr/X11R6/lib 


all: $(OBJ_NAME1) $(BIN_NAME1)
1:   $(OBJ_NAME1) $(BIN_NAME1)


clean:	
	rm -f *.o
	rm -f *.a
	rm -f $(OBJ_NAME1) $(BIN_NAME1)

1c:
	rm -f $(OBJ_NAME1) $(BIN_NAME1)

yos:
	./frontend yos1.pgm yos2.pgm 2 yos_truth.F

ett:
	./frontend ett1.pgm ett2.pgm 2


$(OBJ_NAME1): $(SRC_NAME1) 
	$(CC) $(NORMAL_CFLAGS) $(ARCH_CFLAGS) $(EXTRA_CFLAGS) $(SPECIAL_CFLAGS) -c $(SRC_NAME1) $(LIBS)
$(BIN_NAME1): $(OBJ_NAME1)
	$(CC) $(NORMAL_CFLAGS) $(ARCH_CFLAGS) $(EXTRA_CFLAGS) $(SPECIAL_CFLAGS) -o $(BIN_NAME1) $(OBJ_NAME1) $(LIBS) $(OGL_LIBS_PATH) $(OGL_LIBS) 


