all: main

# Create executable from object file
main: main.o
	gcc main.o -o main

# Create object file from source file
main.o: main.c
	gcc -c main.c

# Clean build files
clean:
	rm -f *.o main