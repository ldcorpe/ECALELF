#!/bin/bash
N=`bjobs | grep  $USER | wc -l`
M=`bjobs| grep RUN | wc -l`
O=`bjobs| grep PEND | wc -l`

echo " Jobs: $N "
echo " RUN : $M"
echo " PEND: $O "
bjobs
