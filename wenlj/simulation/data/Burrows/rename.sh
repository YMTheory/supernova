#!/bin/bash


find . -name "spectra_share*.dat" | xargs sed -i '/^[@#]/d'
