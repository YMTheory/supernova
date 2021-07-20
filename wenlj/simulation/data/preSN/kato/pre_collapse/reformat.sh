#!/bin/bash

cat file | awk -F' ' '{print $2" "$6" "$5" "$1}' > lightcurve_all.dat #nux_bar
cat file | awk -F' ' '{print $2" "$6" "$8" "$1}' > lightcurve_all.dat #nux

rename sum_mu_nu all *
rename sum_mu all *

