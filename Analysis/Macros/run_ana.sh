#!/bin/sh

echo $1

root -b<<EOF
.x ana.C("$1");
EOF
