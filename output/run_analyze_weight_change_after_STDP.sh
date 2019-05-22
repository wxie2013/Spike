#!/bin/sh

root -b<<EOF
.x analyze_weight_change_after_STDP.C++("$1", "$2");
EOF
