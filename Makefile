CFLAGS = -O3  # -funroll-loops
CXXFLAGS = -std=c++17

all:
	g++ $(CFLAGS) $(CXXFLAGS) parser.cpp ppm.cpp raytracer.cpp raytracer_math.cpp tinyxml2.cpp -o raytracer 

test: test.cpp
	g++ -O3 test.cpp -o test -std=c++17 -g
    
png:
	for foo in *.ppm; do convert "$$foo" "$${foo%.ppm}.png"; done
	
