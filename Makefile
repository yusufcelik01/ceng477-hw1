CFLAGS =  -Wall -g  -O3# -funroll-loops
CXXFLAGS = -std=c++17

define RENDER_ALL =
	echo "simple.xml            " >> time.txt ;  (time ./raytracer inputs/simple.xml             ) &>> time.txt; echo "-----------------" >> time.txt;
endef

all:
	g++ $(CFLAGS) $(CXXFLAGS) parser.cpp ppm.cpp raytracer.cpp raytracer_math.cpp tinyxml2.cpp -o raytracer 

test: test.cpp
	g++ -O3 raytracer_math.cpp test.cpp -o test -std=c++17 -g
    
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
	for foo in inputs/*xml; do ./raytracer "$$foo" ;   done

my_important_task: ; $(RENDER_ALL)

.ONESHELL:
