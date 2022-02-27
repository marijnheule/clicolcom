TIMEFORMAT=%R

time ./new-loop.sh graphs/1-FullIns_3.col
time ./new-loop.sh graphs/1-FullIns_4.col
time ./new-loop.sh graphs/1-FullIns_5.col
time ./new-loop.sh graphs/1-Insertions_4.col
#./new-loop.sh graphs/1-Insertions_5.col
#./new-loop.sh graphs/1-Insertions_6.col
time ./new-loop.sh graphs/2-FullIns_3.col
time ./new-loop.sh graphs/2-FullIns_4.col
time ./new-loop.sh graphs/2-FullIns_5.col
time ./new-loop.sh graphs/2-Insertions_3.col
#./new-loop.sh graphs/2-Insertions_4.col
#./new-loop.sh graphs/2-Insertions_5.col
time ./new-loop.sh graphs/3-FullIns_3.col
time ./new-loop.sh graphs/3-FullIns_4.col
time ./new-loop.sh graphs/3-FullIns_5.col
time ./new-loop.sh graphs/3-Insertions_3.col
#./new-loop.sh graphs/3-Insertions_4.col
#./new-loop.sh graphs/3-Insertions_5.col
time ./new-loop.sh graphs/4-FullIns_3.col
time ./new-loop.sh graphs/4-FullIns_4.col
time ./new-loop.sh graphs/4-FullIns_5.col
time ./new-loop.sh graphs/4-Insertions_3.col
#./new-loop.sh graphs/4-Insertions_4.col
time ./new-loop.sh graphs/5-FullIns_3.col
time ./new-loop.sh graphs/5-FullIns_4.col
#./new-loop.sh graphs/C2000.5.col
#./new-loop.sh graphs/C2000.9.col
#./new-loop.sh graphs/C4000.5.col
#./new-loop.sh graphs/DSJC1000.1.col
#./new-loop.sh graphs/DSJC1000.5.col
#./new-loop.sh graphs/DSJC1000.9.col
time ./new-loop.sh graphs/DSJC125.1.col
#./new-loop.sh graphs/DSJC125.5.col
####./new-loop.sh graphs/DSJC125.9.col # seems too difficult
#./new-loop.sh graphs/DSJC250.1.col
#./new-loop.sh graphs/DSJC250.5.col
#./new-loop.sh graphs/DSJC250.9.col
#./new-loop.sh graphs/DSJC500.1.col
#./new-loop.sh graphs/DSJC500.5.col
#./new-loop.sh graphs/DSJC500.9.col
time ./new-loop.sh graphs/DSJR500.1.col
time ./new-loop.sh graphs/DSJR500.1c.col  ## doable but it takes a while
###./new-loop.sh graphs/DSJR500.5.col   ## might be doable. yalsat with --eager gets down to 6 falsified clauses
#time ./new-loop.sh graphs/abb313GPIA.col # depends on the seed
time ./new-loop.sh graphs/anna.col
time ./new-loop.sh graphs/ash331GPIA.col
time ./new-loop.sh graphs/ash608GPIA.col
time ./new-loop.sh graphs/ash958GPIA.col
time ./new-loop.sh graphs/david.col
#./new-loop.sh graphs/flat1000_50_0.col
#./new-loop.sh graphs/flat1000_60_0.col
#./new-loop.sh graphs/flat1000_76_0.col
#./new-loop.sh graphs/flat300_20_0.col
#./new-loop.sh graphs/flat300_26_0.col
#./new-loop.sh graphs/flat300_28_0.col
time ./new-loop.sh graphs/fpsol2.i.1.col
time ./new-loop.sh graphs/fpsol2.i.2.col
time ./new-loop.sh graphs/fpsol2.i.3.col
time ./new-loop.sh graphs/games120.col
time ./new-loop.sh graphs/homer.col
time ./new-loop.sh graphs/huck.col
time ./new-loop.sh graphs/inithx.i.1.col
time ./new-loop.sh graphs/inithx.i.2.col
time ./new-loop.sh graphs/inithx.i.3.col
time ./new-loop.sh graphs/jean.col
#./new-loop.sh graphs/latin_square_10.col
time ./new-loop.sh graphs/le450_15a.col
time ./new-loop.sh graphs/le450_15b.col
time ./new-loop.sh graphs/le450_15c.col  # easy with SLS
time ./new-loop.sh graphs/le450_15d.col  # easy with SLS
time ./new-loop.sh graphs/le450_25a.col
time ./new-loop.sh graphs/le450_25b.col
###./new-loop.sh graphs/le450_25c.col  # possibly with SLS
###./new-loop.sh graphs/le450_25d.col  # possibly with SLS
time ./new-loop.sh graphs/le450_5a.col
time ./new-loop.sh graphs/le450_5b.col
time ./new-loop.sh graphs/le450_5c.col
time ./new-loop.sh graphs/le450_5d.col
time ./new-loop.sh graphs/miles1000.col
time ./new-loop.sh graphs/miles1500.col
time ./new-loop.sh graphs/miles250.col
time ./new-loop.sh graphs/miles500.col
time ./new-loop.sh graphs/miles750.col
time ./new-loop.sh graphs/mug100_1.col
time ./new-loop.sh graphs/mug100_25.col
time ./new-loop.sh graphs/mug88_1.col
time ./new-loop.sh graphs/mug88_25.col
time ./new-loop.sh graphs/mulsol.i.1.col
time ./new-loop.sh graphs/mulsol.i.2.col
time ./new-loop.sh graphs/mulsol.i.3.col
time ./new-loop.sh graphs/mulsol.i.4.col
time ./new-loop.sh graphs/mulsol.i.5.col
time ./new-loop.sh graphs/myciel3.col
time ./new-loop.sh graphs/myciel4.col
time ./new-loop.sh graphs/myciel5.col
time ./new-loop.sh graphs/myciel6.col
#./new-loop.sh graphs/myciel7.col
#time ./new-loop.sh graphs/qg.order100.col #doable but expensive
time ./new-loop.sh graphs/qg.order30.col
time ./new-loop.sh graphs/qg.order40.col #doable but expensive
time ./new-loop.sh graphs/qg.order60.col
#./new-loop.sh graphs/queen10_10.col
#./new-loop.sh graphs/queen11_11.col
#./new-loop.sh graphs/queen12_12.col
#./new-loop.sh graphs/queen13_13.col
#./new-loop.sh graphs/queen14_14.col
#./new-loop.sh graphs/queen15_15.col
#./new-loop.sh graphs/queen16_16.col
time ./new-loop.sh graphs/queen5_5.col
time ./new-loop.sh graphs/queen6_6.col
time ./new-loop.sh graphs/queen7_7.col
time ./new-loop.sh graphs/queen8_12.col
time ./new-loop.sh graphs/queen8_8.col
#time ./new-loop.sh graphs/queen9_9.col  # can be hard
time ./new-loop.sh graphs/r1000.1.col
#./new-loop.sh graphs/r1000.1c.col
#./new-loop.sh graphs/r1000.5.col
time ./new-loop.sh graphs/r125.1.col
time ./new-loop.sh graphs/r125.1c.col
time ./new-loop.sh graphs/r125.5.col
time ./new-loop.sh graphs/r250.1.col
time ./new-loop.sh graphs/r250.1c.col
time ./new-loop.sh graphs/r250.5.col
time ./new-loop.sh graphs/school1.col
time ./new-loop.sh graphs/school1_nsh.col
time ./new-loop.sh graphs/wap01a.col
time ./new-loop.sh graphs/wap02a.col
#./new-loop.sh graphs/wap03a.col
#./new-loop.sh graphs/wap04a.col
time ./new-loop.sh graphs/wap05a.col
time ./new-loop.sh graphs/wap06a.col
#./new-loop.sh graphs/wap07a.col
#./new-loop.sh graphs/wap08a.col
time ./new-loop.sh graphs/will199GPIA.col
time ./new-loop.sh graphs/zeroin.i.1.col
time ./new-loop.sh graphs/zeroin.i.2.col
time ./new-loop.sh graphs/zeroin.i.3.col
