all:
	g++ -O3 parser.cpp ppm.cpp raytracer.cpp tinyxml2.cpp -o raytracer -std=c++11 -g

test: test.cpp
	g++ -O3 test.cpp -o test -std=c++11 -g
    
png:
	for foo in *.ppm; do convert "$$foo" "$${foo%.ppm}.png"; done
	
