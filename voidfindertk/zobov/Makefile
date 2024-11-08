
.PHONY : loader

loader:
	gcc -fPIC -shared -o zobov_loader.so zobov_loader.c
	gcc -fPIC -shared -o zones_in_void.so zones_in_void.c
	gcc -fPIC -shared -o tracers_in_zones.so tracers_in_zones.c


clean:
	make -C clean

