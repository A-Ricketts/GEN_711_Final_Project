#! /bin/bash

source activate genomics

data=$1
altdata=$(ls $data)

ls $data | grep -vo "_R*" | echo
