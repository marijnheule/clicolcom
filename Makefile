all: color filter maxclique

color: color.c
	gcc color.c -O2 -o color

filter: filter.c
	gcc filter.c -O2 -o filter

maxclique: maxclique.c
	gcc maxclique.c -O2 -o maxclique

clean:
	rm color filter maxclique
