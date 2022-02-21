all: color filter maxclique optimize-ord inc-max-clique

color: color.c
	gcc color.c -O2 -o color

filter: filter.c
	gcc filter.c -O2 -o filter

maxclique: maxclique.c
	gcc maxclique.c -O2 -o maxclique

optimize-ord: optimize-ord.c
	gcc optimize-ord.c -O2 -o optimize-ord

inc-max-clique: NewIncMaxCLQictai13.c
	gcc -Wall NewIncMaxCLQictai13.c -O2 -o inc_max_clique

clean:
	rm color filter maxclique optimize-ord inc_max_clique
