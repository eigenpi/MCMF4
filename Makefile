CC = g++

CLASSDIR = /home/cristinel/projects/mcmf4
INCDIRS = $(CLASSDIR)/
LIB_DIR = -L/usr/lib/X11
LIB = -lX11 -lm

#OPT_FLAGS =
OPT_FLAGS = -O3
WARN_FLAGS =
#WARN_FLAGS = -Wall -Wpointer-arith -Wcast-qual -Wstrict-prototypes -O -D__USE_FIXED_PROTOTYPES__ -ansi -pedantic -Wmissing-prototypes -Wshadow -Wcast-align -D_POSIX_SOURCE
DEBUG_FLAGS =
#DEBUG_FLAGS = -g

FLAGS = $(OPT_FLAGS)
FLAGS += $(WARN_FLAGS)
FLAGS += $(DEBUG_FLAGS)
FLAGS += $(addprefix -I, $(INCDIRS))

EXE = mcmf4
OBJ = mcmf4.o 
SRC = mcmf4.cpp
H = mcmf4.h 

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) -o $(EXE) $(LIB_DIR) $(LIB)

mcmf4.o: mcmf4.cpp $(H)
	$(CC) -c $(FLAGS) mcmf4.cpp
