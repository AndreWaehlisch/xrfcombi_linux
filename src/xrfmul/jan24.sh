#! /bin/bash
cp ./data/S5min.lin ./data/lines.dat
cp ./data/S5min.dat ./data/sample.dat
xrfmul > res.S5minLa2
cp ./data/SYSZ.lin ./data/lines.dat
cp ./data/SYSZ.dat ./data/sample.dat
xrfmul > res.SYSZLa2
cp ./data/S5mindup.lin ./data/lines.dat
cp ./data/S5mindup.dat ./data/sample.dat
xrfmul > res.S5mindupLa2
echo "Ready with script!!"

