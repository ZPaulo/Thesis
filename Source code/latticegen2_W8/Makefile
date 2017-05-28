CC=g++
CFLAGS=
LDFLAGS=
LIBS=-lm

CPP_FILES=main.cpp FluidTracking.cpp GeomRead.cpp Intersect.cpp MeshGen.cpp Vector2D.cpp Vector3D.cpp
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CPP_FILES:.cpp=.o))
OBJ_DIR=$(PRE_NAME)obj

EXEC=latticegen

all: $(EXEC)
$(EXEC): $(OBJ_DIR) $(OBJ_FILES); \
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $(OBJ_FILES)

$(OBJ_DIR):; \
	mkdir $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.cpp; \
	$(CC) $(CFLAGS) $(LIBS) -c $< -o $@

.PHONY: clean

clean:; \
	rm -rf $(OBJ_DIR) $(EXEC)
