CC = g++
CFLAGS = -g -Wall -O2 -march=znver1 -pipe -std=c++1z 
target = lock_in

all: $(target).cpp 
	$(CC) $(CFLAGS) -o $(target) $(target).cpp -fopenmp
	./$(target)

clean: 
	$(RM) myprog