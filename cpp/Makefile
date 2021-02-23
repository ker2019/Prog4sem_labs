GMSH_LIB_PATH=/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/lib
GMSH_INCLUDE_PATH=/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/include

.PHONY: all clean

all: torus heart

torus: torus.cpp
	g++ -o torus torus.cpp -lgmsh -L${GMSH_LIB_PATH} -I${GMSH_INCLUDE_PATH}

heart: heart.cpp
	g++ -o heart heart.cpp -lgmsh -L${GMSH_LIB_PATH} -I${GMSH_INCLUDE_PATH}

clean:
	rm -f torus heart
