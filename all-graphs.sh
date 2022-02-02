TIMEFORMAT=%R

time ./loop-solve.sh graphs/1-FullIns_3.col
time ./loop-solve.sh graphs/1-FullIns_4.col
time ./loop-solve.sh graphs/1-FullIns_5.col
time ./loop-solve.sh graphs/1-Insertions_4.col
#./loop-solve.sh graphs/1-Insertions_5.col
#./loop-solve.sh graphs/1-Insertions_6.col
time ./loop-solve.sh graphs/2-FullIns_3.col
time ./loop-solve.sh graphs/2-FullIns_4.col
time ./loop-solve.sh graphs/2-FullIns_5.col
time ./loop-solve.sh graphs/2-Insertions_3.col
#./loop-solve.sh graphs/2-Insertions_4.col
#./loop-solve.sh graphs/2-Insertions_5.col
time ./loop-solve.sh graphs/3-FullIns_3.col
time ./loop-solve.sh graphs/3-FullIns_4.col
time ./loop-solve.sh graphs/3-FullIns_5.col
time ./loop-solve.sh graphs/3-Insertions_3.col
#./loop-solve.sh graphs/3-Insertions_4.col
#./loop-solve.sh graphs/3-Insertions_5.col
time ./loop-solve.sh graphs/4-FullIns_3.col
time ./loop-solve.sh graphs/4-FullIns_4.col
#./loop-solve.sh graphs/4-FullIns_5.col
time ./loop-solve.sh graphs/4-Insertions_3.col
#./loop-solve.sh graphs/4-Insertions_4.col
time ./loop-solve.sh graphs/5-FullIns_3.col
time ./loop-solve.sh graphs/5-FullIns_4.col
#./loop-solve.sh graphs/C2000.5.col
#./loop-solve.sh graphs/C2000.9.col
#./loop-solve.sh graphs/C4000.5.col
#./loop-solve.sh graphs/DSJC1000.1.col
#./loop-solve.sh graphs/DSJC1000.5.col
#./loop-solve.sh graphs/DSJC1000.9.col
time ./loop-solve.sh graphs/DSJC125.1.col
#./loop-solve.sh graphs/DSJC125.5.col
####./loop-solve.sh graphs/DSJC125.9.col # seems too difficult
#./loop-solve.sh graphs/DSJC250.1.col
#./loop-solve.sh graphs/DSJC250.5.col
#./loop-solve.sh graphs/DSJC250.9.col
#./loop-solve.sh graphs/DSJC500.1.col
#./loop-solve.sh graphs/DSJC500.5.col
#./loop-solve.sh graphs/DSJC500.9.col
time ./loop-solve.sh graphs/DSJR500.1.col
#time ./loop-solve.sh graphs/DSJR500.1c.col  ## doable but it takes a while
###./loop-solve.sh graphs/DSJR500.5.col   ## might be doable. yalsat with --eager gets down to 6 falsified clauses
#time ./loop-solve.sh graphs/abb313GPIA.col # depends on the seed
time ./loop-solve.sh graphs/anna.col
time ./loop-solve.sh graphs/ash331GPIA.col
time ./loop-solve.sh graphs/ash608GPIA.col
time ./loop-solve.sh graphs/ash958GPIA.col
time ./loop-solve.sh graphs/david.col
#./loop-solve.sh graphs/flat1000_50_0.col
#./loop-solve.sh graphs/flat1000_60_0.col
#./loop-solve.sh graphs/flat1000_76_0.col
#./loop-solve.sh graphs/flat300_20_0.col
#./loop-solve.sh graphs/flat300_26_0.col
#./loop-solve.sh graphs/flat300_28_0.col
time ./loop-solve.sh graphs/fpsol2.i.1.col
time ./loop-solve.sh graphs/fpsol2.i.2.col
time ./loop-solve.sh graphs/fpsol2.i.3.col
time ./loop-solve.sh graphs/games120.col
time ./loop-solve.sh graphs/homer.col
time ./loop-solve.sh graphs/huck.col
time ./loop-solve.sh graphs/inithx.i.1.col
time ./loop-solve.sh graphs/inithx.i.2.col
time ./loop-solve.sh graphs/inithx.i.3.col
time ./loop-solve.sh graphs/jean.col
#./loop-solve.sh graphs/latin_square_10.col
time ./loop-solve.sh graphs/le450_15a.col
time ./loop-solve.sh graphs/le450_15b.col
time ./loop-solve.sh graphs/le450_15c.col  # easy with SLS
time ./loop-solve.sh graphs/le450_15d.col  # easy with SLS
time ./loop-solve.sh graphs/le450_25a.col
time ./loop-solve.sh graphs/le450_25b.col
###./loop-solve.sh graphs/le450_25c.col  # possibly with SLS
###./loop-solve.sh graphs/le450_25d.col  # possibly with SLS
time ./loop-solve.sh graphs/le450_5a.col
time ./loop-solve.sh graphs/le450_5b.col
time ./loop-solve.sh graphs/le450_5c.col
time ./loop-solve.sh graphs/le450_5d.col
time ./loop-solve.sh graphs/miles1000.col
time ./loop-solve.sh graphs/miles1500.col
time ./loop-solve.sh graphs/miles250.col
time ./loop-solve.sh graphs/miles500.col
time ./loop-solve.sh graphs/miles750.col
time ./loop-solve.sh graphs/mug100_1.col
time ./loop-solve.sh graphs/mug100_25.col
time ./loop-solve.sh graphs/mug88_1.col
time ./loop-solve.sh graphs/mug88_25.col
time ./loop-solve.sh graphs/mulsol.i.1.col
time ./loop-solve.sh graphs/mulsol.i.2.col
time ./loop-solve.sh graphs/mulsol.i.3.col
time ./loop-solve.sh graphs/mulsol.i.4.col
time ./loop-solve.sh graphs/mulsol.i.5.col
time ./loop-solve.sh graphs/myciel3.col
time ./loop-solve.sh graphs/myciel4.col
time ./loop-solve.sh graphs/myciel5.col
time ./loop-solve.sh graphs/myciel6.col
#./loop-solve.sh graphs/myciel7.col
#time ./loop-solve.sh graphs/qg.order100.col #doable but expensive
time ./loop-solve.sh graphs/qg.order30.col
time ./loop-solve.sh graphs/qg.order40.col #doable but expensive
#time ./loop-solve.sh graphs/qg.order60.col
#./loop-solve.sh graphs/queen10_10.col
#./loop-solve.sh graphs/queen11_11.col
#./loop-solve.sh graphs/queen12_12.col
#./loop-solve.sh graphs/queen13_13.col
#./loop-solve.sh graphs/queen14_14.col
#./loop-solve.sh graphs/queen15_15.col
#./loop-solve.sh graphs/queen16_16.col
time ./loop-solve.sh graphs/queen5_5.col
time ./loop-solve.sh graphs/queen6_6.col
time ./loop-solve.sh graphs/queen7_7.col
time ./loop-solve.sh graphs/queen8_12.col
time ./loop-solve.sh graphs/queen8_8.col
#time ./loop-solve.sh graphs/queen9_9.col  # can be hard
time ./loop-solve.sh graphs/r1000.1.col
#./loop-solve.sh graphs/r1000.1c.col
#./loop-solve.sh graphs/r1000.5.col
time ./loop-solve.sh graphs/r125.1.col
time ./loop-solve.sh graphs/r125.1c.col
time ./loop-solve.sh graphs/r125.5.col
time ./loop-solve.sh graphs/r250.1.col
time ./loop-solve.sh graphs/r250.1c.col
time ./loop-solve.sh graphs/r250.5.col
time ./loop-solve.sh graphs/school1.col
time ./loop-solve.sh graphs/school1_nsh.col
time ./loop-solve.sh graphs/wap01a.col
time ./loop-solve.sh graphs/wap02a.col
#./loop-solve.sh graphs/wap03a.col
#./loop-solve.sh graphs/wap04a.col
time ./loop-solve.sh graphs/wap05a.col
time ./loop-solve.sh graphs/wap06a.col
#./loop-solve.sh graphs/wap07a.col
#./loop-solve.sh graphs/wap08a.col
time ./loop-solve.sh graphs/will199GPIA.col
time ./loop-solve.sh graphs/zeroin.i.1.col
time ./loop-solve.sh graphs/zeroin.i.2.col
time ./loop-solve.sh graphs/zeroin.i.3.col
