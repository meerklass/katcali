#! /bin/bash
echo "# data list creating"

ls *h_*data >h.dat
echo "HH lines: "
wc h.dat
echo "VV lines: "
ls *v_*data >v.dat
wc v.dat

echo "# list reshape and compare"
echo "HH: "
awk -F _ '{print $1,$2}' h.dat | sed -e 's/h$//g' > result_level3_hlist.txt
wc result_level3_hlist.txt

echo "--------------------------------"
tail -10 h.dat
echo "--------------------------------"
tail -10 result_level3_hlist.txt
echo "--------------------------------"
echo "--------------------------------"

echo "VV: "
awk -F _ '{print $1,$2}' v.dat | sed -e 's/v$//g' > result_level3_vlist.txt
wc result_level3_vlist.txt

echo "--------------------------------"
tail -10 v.dat
echo "--------------------------------"
tail -10 result_level3_vlist.txt

echo "#"


