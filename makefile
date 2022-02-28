CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl -lm -lgslcblas
#---------------------------------------------------------
$(OBJ)/sim1.o:\
		sim1.C
	$(CC) $(CFLAGS) -o $(OBJ)/sim1.o -c \
		sim1.C
#---------------------------------------------------------
sim1: \
		$(OBJ)/sim1.o
	$(CC) $(LDFLAGS) -o sim1 \
		$(OBJ)/sim1.o \
		$(LIBS)
#---------------------------------------------
