#!/bin/bash

file=$argv$1

awk '{print $1,$2}' ${file} 
