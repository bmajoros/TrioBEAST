CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -I$(GSLDIR)/include -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl -lm -lgslcblas
#---------------------------------------------------------
$(OBJ)/sim2.o:\
		sim2.C
	$(CC) $(CFLAGS) -o $(OBJ)/sim2.o -c \
		sim2.C
#---------------------------------------------------------
sim2: \
		$(OBJ)/sim2.o
	$(CC) $(LDFLAGS) -o sim2 \
		$(OBJ)/sim2.o \
		$(LIBS)
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
#--------------------------------------------------------
$(OBJ)/phase-trio.o:\
		phase-trio.C
	$(CC) $(CFLAGS) -o $(OBJ)/phase-trio.o -c \
		phase-trio.C
#---------------------------------------------------------
phase-trio: \
		$(OBJ)/phase-trio.o
	$(CC) $(LDFLAGS) -o phase-trio \
		$(OBJ)/phase-trio.o \
		$(LIBS)
#---------------------------------------------------------
