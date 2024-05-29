
all: qhull voboz

.PHONY : qhull voboz examples
.PHONY : loader

loader:
	gcc -fPIC -shared -o zobov_loader.so zobov_loader.c
	g++ -o lectura_zones lectura_zones_in_void.cpp
	g++ -o particle_zones lectura_zonfile.cpp

qhull:
	make -C qhull/src

voboz: 
	make -C src install

examples:
	make -C examples
	make -C examples data

clean:
	make -C src clean
	make -C qhull/src cleanall
	#make -C examples clean
	rm -f bin/*
	rm -f zobov_loader.so
