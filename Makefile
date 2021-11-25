CFLAGS =  -Wall -g  -O3# -funroll-loops
CXXFLAGS = -std=c++17

all:
	g++ $(CFLAGS) $(CXXFLAGS) parser.cpp ppm.cpp raytracer.cpp raytracer_math.cpp tinyxml2.cpp -o raytracer 

test: test.cpp
	g++ -O3 raytracer_math.cpp test.cpp -o test -std=c++17 -g
    
png:
	for foo in *.ppm; do convert "$$foo" "$${foo%.ppm}.png"; done

bash:
	for foo in inputs/*.xml; do ./bash_test.sh ; echo "input/$$foo"; done
	
#render_small:
#	-./raytracer inputs/simple.xml
#	-./raytracer inputs/simple_shading.xml
#	-./raytracer inputs/simple_reflectance.xml
#	-./raytracer inputs/mirror_spheres.xml

#render_small:
get_time:
	echo "simple.xml            " >> time.txt ;  (time ./raytracer inputs/simple.xml             ) &>> time.txt; echo "-----------------" >> time.txt;
#echo "simple_shading.xml    " >> time.txt ;  (time ./raytracer inputs/simple_shading.xml     ) &>> time.txt; echo "-----------------" >> time.txt;
#echo "simple_reflectance.xml" >> time.txt ;  (time ./raytracer inputs/simple_reflectance.xml ) &>> time.txt; echo "-----------------" >> time.txt;
#echo "mirror_spheres.xml    " >> time.txt ;  (time ./raytracer inputs/mirror_spheres.xml     ) &>> time.txt; echo "-----------------" >> time.txt;


render__all:
	for foo in inputs/*xml; echo "$$foo" >> time.txt ; do (time ./raytracer "$$foo") &>> time.txt; echo "-----------------" >> time.txt;   done

