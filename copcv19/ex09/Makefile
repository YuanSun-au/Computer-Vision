CC = g++

NAME1 = frontend
NAME2 = view_epi

BIN_NAME1 = $(NAME1)
OBJ_NAME1 = $(NAME1).o
SRC_NAME1 = $(NAME1).cc

BIN_NAME2 = $(NAME2)
OBJ_NAME2 = $(NAME2).o
SRC_NAME2 = $(NAME2).cc

NORMAL_CFLAGS = -g -funroll-loops -w
#EXTRA_CFLAGS = -fomit-frame-pointer #-static


INCL_PATH =

LIBS_PATH =
LIBS = #-lm

OGL_LIBS = -lGL -lGLU -lglut 
# -lXxf86vm 
OGL_LIBS_PATH = -L/usr/X11R6/lib 


all: $(OBJ_NAME1) $(BIN_NAME1) $(OBJ_NAME2) $(BIN_NAME2)
1:   $(OBJ_NAME1) $(BIN_NAME1)
2:   $(OBJ_NAME2) $(BIN_NAME2)


runtsu:
	./frontend tsu1.pgm tsu2.pgm 1
viewtsu:
	./view_epi tsu1.pgm tsu2.pgm tsu_t.fm matrix_e.fm 1

runjav:
	./frontend jav1.pgm jav2.pgm 1
viewjav:
	./view_epi jav1.pgm jav2.pgm jav_t.fm matrix_e.fm 1


clean:	
	rm -f *.o
	rm -f *.a
	rm -f $(OBJ_NAME1) $(BIN_NAME1)
	rm -f $(OBJ_NAME2) $(BIN_NAME2)	

1c:
	rm -f $(OBJ_NAME1) $(BIN_NAME1)
2c:
	rm -f $(OBJ_NAME2) $(BIN_NAME2)


$(OBJ_NAME1): $(SRC_NAME1) 
	$(CC) $(NORMAL_CFLAGS) $(ARCH_CFLAGS) $(EXTRA_CFLAGS) $(SPECIAL_CFLAGS) -c $(SRC_NAME1) $(LIBS)
$(BIN_NAME1): $(OBJ_NAME1)
	$(CC) $(NORMAL_CFLAGS) $(ARCH_CFLAGS) $(EXTRA_CFLAGS) $(SPECIAL_CFLAGS) -o $(BIN_NAME1) $(OBJ_NAME1) $(LIBS) $(OGL_LIBS_PATH) $(OGL_LIBS) 
$(OBJ_NAME2): $(SRC_NAME2) 
	$(CC) $(NORMAL_CFLAGS) $(ARCH_CFLAGS) $(EXTRA_CFLAGS) $(SPECIAL_CFLAGS) -c $(SRC_NAME2) $(LIBS)
$(BIN_NAME2): $(OBJ_NAME2)
	$(CC) $(NORMAL_CFLAGS) $(ARCH_CFLAGS) $(EXTRA_CFLAGS) $(SPECIAL_CFLAGS) -o $(BIN_NAME2) $(OBJ_NAME2) $(LIBS) $(OGL_LIBS_PATH) $(OGL_LIBS) 


