SHELL := /bin/bash
CXX=g++
CFLAGS =  -g  -O3  # -funroll-loops
CXXFLAGS = -std=c++17
OBJECT_FILES= parser.o ppm.o raytracer_math.cpp tinyxml2.o 


all: $(OBJECT_FILES)
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(OBJECT_FILES)  notth_raytracer.cpp  -o raytracer 



%.o: %.cpp %.h
	$(CXX)  $(CFLAGS) $(CXXFLAGS) $< -c

thread:
	$(CXX) $(CFLAGS) $(CXXFLAGS) parser.cpp ppm.cpp raytracer.cpp raytracer_math.cpp tinyxml2.cpp -o raytracer  -pthread

test: test.cpp
	$(CXX) -O3 raytracer_math.cpp test.cpp -o test -std=c++17 -g
    
png:
	for foo in *.ppm; do convert "$$foo" "$${foo%.ppm}.png"; done

bash:
	for foo in inputs/*.xml; do ./bash_test.sh ; echo "input/$$foo"; done
	
render_small:
	-./raytracer inputs/simple.xml
	-./raytracer inputs/simple_shading.xml
	-./raytracer inputs/simple_reflectance.xml
	-./raytracer inputs/mirror_spheres.xml



render__all:
	for foo in inputs/*xml; do bar="time_$${foo#inputs/}"; (time ./raytracer "$$foo") 2> "$${bar%.xml}.txt" ;  done


