CC = g++

_EXE = main
_OBJ = frontend.o
_SRC = frontend.cc
_INCL = 

# correct folders
ODIR = ../obj
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
SRCDIR = .
SRC = $(patsubst %,$(SRCDIR)/%,$(_SRC))
INCLDIR = ../include
INCL = $(patsubst %,$(INCDIR)/%,$(_INCL))
EXEDIR = ..
EXE = $(patsubst %,$(EXEDIR)/%,$(_EXE))
LIBDIR = ../lib


CFLAGS = -funroll-loops -w -I$(INCLDIR) -Wall
OPTIMIZATION = -O3
EXTRA_CFLAGS = #-fomit-frame-pointer #-static

# libraries
LIBS = -L$(LIBDIR) #-lm
OGL_LIBS = -lGL -lGLU -lglut #-lXxf86vm 
OGL_LIBS_PATH = -L/usr/X11R6/lib 



$(EXE): $(OBJ)
	$(CC) $^ -o $@ $(LIBS) $(OGL_LIBS_PATH) $(OGL_LIBS) $(CFLAGS) $(OPTIMIZATION) $(EXTRA_CFLAGS) $(DBGCFLAGS) 

# $(OBJ): $(SRC) $(INCL)
$(ODIR)/%.o: $(SRCDIR)/%.cc $(INCL)
	@mkdir -p $(@D)
	$(CC) $< -o $@ -c $(CFLAGS) $(OPTIMIZATION) $(EXTRA_CFLAGS) $(DBGCFLAGS)



all: $(OBJ) $(EXE)
debug: DBGCFLAGS = -DDEBUG -g
debug: OPTIMIZATION = -O0			# optimization level = 0 (otherwise breakpoints will not work with optimized code)
debug: $(OBJ) $(EXE)

clean:	
	rm -f *.o
	rm -f *.a
	rm -f $(OBJ) $(EXE)

1c:
	rm -f $(OBJ) $(EXE)

.PHONY: all clean lc debug