#
# Makefile to build up the several programs of the package. it
# wouldn't be a bad idea to have a look at the file README.1st...
#
#
# in the next lines you have to choose the right parameters for your
# system. I've put as default values the ones corresponding to the GNU
# C/C++ compiler, but you should change them if you want to use a
# different compiler. I've added (in commented lines) the corresponding
# values for the compiler that comes with HP UX 10.20.
#
# Ansi C compiler
#
CC=gcc
# CC=cc (this is for HP UX 10.20)
#
# C++ compiler
#
CP=g++
# CP=CC (this is for HP UX 10.20)
#
# compilation flags for the C compiler
#
CFLAGS=-O3 -Wall
# CFLAGS=-O -Aa (this is for HP UX 10.20)
#
# compilation flags for the C++ compiler
#
CPPFLAGS=-O3 -Wall
# CPPFLAGS=-O (this is for HP UX 10.20)
#
# linking flags
#
LFLAGS=-s -lm
#
# directories to store the binaries and to find/store data files resp.
# I've used relative pathnames, but it is better to use absolute ones
# (with the relative pathnames, you have to execute the programs from
# the BIN directory, otherwise they could not find the DATA directory).
#
# BIN=../bin/
BIN=C:/Users/19525/OneDriveUCB/Documents/Toolkit/bin/
#DAT=../data/
DAT=C:/Users/19525/OneDriveUCB/Documents/Toolkit/data/
#
# directory to put working files (they will be erased at the end of
# the execution).
#
#TMP=../data/
TMP=C:/Users/19525/OneDriveUCB/Documents/Toolkit/data/
#
# ========================================================
#  you shouldn't need to modify anything beyond this line
# ========================================================
#
#######################################################
# to build up the programs that transform coordinates #
#######################################################

TCNF=main-tcnf.o repo.o erpe6.o io4s.o io6s.o io6p.o mp4s.o mp6s.o mp6p.o iol.o litr5.o uneix.o msgs.o

tcnf: $(TCNF)
	$(CC) $(CFLAGS) $(TCNF) -o $(BIN)tcnf $(LFLAGS)
 
#######
# all #
#######

all: tcnf

############################
# how to make the .o files #
############################

main-tcnf.o: main-tcnf.c
	$(CC) -c $(CFLAGS) -DDATA=\"$(DAT)\" main-tcnf.c
erpe6.o: erpe6.c
	$(CC) -c $(CFLAGS) erpe6.c
io6s.o: io6s.c
	$(CC) -c $(CFLAGS) io6s.c
io6p.o: io6p.c
	$(CC) -c $(CFLAGS) io6p.c
io4s.o: io4s.c
	$(CC) -c $(CFLAGS) io4s.c
iol.o: iol.c
	$(CC) -c $(CFLAGS) iol.c
mp6s.o: mp6s.c
	$(CC) -c $(CFLAGS) mp6s.c
mp6p.o: mp6p.c
	$(CC) -c $(CFLAGS) mp6p.c
mp4s.o: mp4s.c
	$(CC) -c $(CFLAGS) mp4s.c
uneix.o: uneix.c
	$(CC) -c $(CFLAGS) uneix.c
msgs.o: msgs.c
	$(CC) -c $(CFLAGS) msgs.c
repo.o: repo.c
	$(CC) -c $(CFLAGS) repo.c
litr5.o: litr5.c
	$(CC) -c $(CFLAGS) litr5.c
#########
# clean #
#########
clean:
#	rm -f *.o
#	rm -f $(BIN)*
	-del *.o
	-del $(BIN)*