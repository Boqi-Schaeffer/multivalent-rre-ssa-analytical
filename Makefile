CC=g++
OPTS= -Wall --pedantic -std=c++17
# LIBS= -lsfml-graphics -lsfml-window -lsfml-system
SRC=SpatialSSA.cpp

start: $(SRC:.cpp=.o)
	$(CC) $(OPTS) -o SpatialSSA $(SRC:.cpp=.o)
circle.o: circle.cpp

# $(CC) $(OPTS) -o SpatialSSA $(SRC:.cpp=.o) $(LIBS)

%.o: %.cpp
	$(CC) $(OPTS) -c $< -o $@

clean:
	rm -f *.o
	rm -f SpatialSSA.exe
