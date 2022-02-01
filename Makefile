all: color filter maxclique optimize-ord

color: color.c
	gcc color.c -O2 -o color

filter: filter.c
	gcc filter.c -O2 -o filter

maxclique: maxclique.c
	gcc maxclique.c -O2 -o maxclique

optimize-ord: optimize-ord.c
	gcc optimize-ord.c -O2 -o optimize-ord

clean:
	rm color filter maxclique optimize-ord
