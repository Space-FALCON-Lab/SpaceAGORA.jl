# Revision 24 target orbits (circular; radius = closing apoapsis radius)

Aerobraking and plane costs DERIVED from the Revision 20 DP (S-fixed, free arrival, node free); preview tables below the Odyssey altitude; V0 at r_p = 3,504 km (ASSUMED). Provenance per entry in target_orbits.json.

| target | i (deg) | alt (km) | r (km) | e final | circle match V0 (m/s) | best geometry (S-fixed) | aero (m/s) | plane (m/s) | total = aero + plane + V0 | Odyssey geometry total | class |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ODY_500km | 93.522 | 500 | 3896 | 0.053 | 89.1 | G1038 (i 93.5) | 0.187 | 0.009 | 89.2 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| ODY_1000km | 93.522 | 1000 | 4396 | 0.113 | 181.5 | G1128 (i 93.5) | 0.129 | 0.004 | 181.7 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| ODY_2000km | 93.522 | 2000 | 5396 | 0.213 | 317.3 | G0842 (i 93.5) | 0.102 | 0.002 | 317.5 | 319.8 | cheap (aerobraking + plane <= 0.5 m/s) |
| ODY_3000km | 93.522 | 3000 | 6396 | 0.292 | 410.5 | G0844 (i 93.5) | 0.080 | 0.001 | 410.6 | 412.0 | cheap (aerobraking + plane <= 0.5 m/s) |
| ODY_ODY | 93.522 | 4012 | 7408 | 0.358 | 477.5 | G0942 (i 93.5) | 0.079 | 0.001 | 477.6 | 478.7 | cheap (aerobraking + plane <= 0.5 m/s) |
| HIGH_500km | 80.000 | 500 | 3896 | 0.053 | 89.1 | G3436 (i 80.0) | 0.071 | 0.000 | 89.1 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| HIGH_1000km | 80.000 | 1000 | 4396 | 0.113 | 181.5 | G3436 (i 80.0) | 0.071 | 0.000 | 181.6 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| HIGH_2000km | 80.000 | 2000 | 5396 | 0.213 | 317.3 | G3432 (i 80.0) | 0.009 | 0.000 | 317.4 | 703.2 | cheap (aerobraking + plane <= 0.5 m/s) |
| HIGH_3000km | 80.000 | 3000 | 6396 | 0.292 | 410.5 | G3432 (i 80.0) | 0.000 | 0.000 | 410.5 | 695.0 | cheap (aerobraking + plane <= 0.5 m/s) |
| HIGH_ODY | 80.000 | 4012 | 7408 | 0.358 | 477.5 | G3366 (i 80.0) | 0.001 | 0.000 | 477.5 | 697.5 | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT63_500km | 63.400 | 500 | 3896 | 0.053 | 89.1 | G0486 (i 63.4) | 0.000 | 0.000 | 89.1 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT63_1000km | 63.400 | 1000 | 4396 | 0.113 | 181.5 | G0486 (i 63.4) | 0.000 | 0.000 | 181.5 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT63_2000km | 63.400 | 2000 | 5396 | 0.213 | 317.3 | G0486 (i 63.4) | 0.000 | 0.000 | 317.3 | 1417.6 | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT63_3000km | 63.400 | 3000 | 6396 | 0.292 | 410.5 | G0486 (i 63.4) | 0.000 | 0.000 | 410.5 | 1301.5 | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT63_ODY | 63.400 | 4012 | 7408 | 0.358 | 477.5 | G0460 (i 63.4) | 0.000 | 0.000 | 477.5 | 1217.5 | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT117_500km | 116.600 | 500 | 3896 | 0.053 | 89.1 | G1230 (i 116.6) | 0.000 | 0.000 | 89.1 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT117_1000km | 116.600 | 1000 | 4396 | 0.113 | 181.5 | G1230 (i 116.6) | 0.000 | 0.000 | 181.5 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT117_2000km | 116.600 | 2000 | 5396 | 0.213 | 317.3 | G1230 (i 116.6) | 0.000 | 0.000 | 317.3 | 1110.5 | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT117_3000km | 116.600 | 3000 | 6396 | 0.292 | 410.5 | G1156 (i 116.6) | 0.000 | 0.000 | 410.5 | 1036.0 | cheap (aerobraking + plane <= 0.5 m/s) |
| CRIT117_ODY | 116.600 | 4012 | 7408 | 0.358 | 477.5 | G1226 (i 116.6) | 0.000 | 0.000 | 477.5 | 985.8 | cheap (aerobraking + plane <= 0.5 m/s) |
| EQ_500km | 20.000 | 500 | 3896 | 0.053 | 89.1 | G0050 (i 20.0) | 0.000 | 0.000 | 89.1 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| EQ_1000km | 20.000 | 1000 | 4396 | 0.113 | 181.5 | G0050 (i 20.0) | 0.000 | 0.000 | 181.5 | nan | cheap (aerobraking + plane <= 0.5 m/s) |
| EQ_2000km | 20.000 | 2000 | 5396 | 0.213 | 317.3 | G0050 (i 20.0) | 0.000 | 0.000 | 317.3 | 3194.7 | cheap (aerobraking + plane <= 0.5 m/s) |
| EQ_3000km | 20.000 | 3000 | 6396 | 0.292 | 410.5 | G0050 (i 20.0) | 0.000 | 0.000 | 410.5 | 2872.0 | cheap (aerobraking + plane <= 0.5 m/s) |
| EQ_ODY | 20.000 | 4012 | 7408 | 0.358 | 477.5 | G0342 (i 20.0) | 0.000 | 0.000 | 477.5 | 2621.5 | cheap (aerobraking + plane <= 0.5 m/s) |
