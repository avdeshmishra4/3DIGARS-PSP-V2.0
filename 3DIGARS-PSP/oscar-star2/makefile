DD = g++ -O
SOURCE = ./source
LIB = ./library
OBJECT = $(LIB)/side4.o $(LIB)/main.o

Main: $(OBJECT)
	$(DD)  $(OBJECT) -lm -o oscar-star

$(LIB)/side4.o: $(SOURCE)/side.h $(SOURCE)/side4.cpp
	$(DD)  $(SOURCE)/side4.cpp -c -o $(LIB)/side4.o

$(LIB)/main.o:  $(SOURCE)/side.h $(SOURCE)/main.cpp
	$(DD)  $(SOURCE)/main.cpp -c -o $(LIB)/main.o
