#!/bin/bash

# create water.xyz
cat > water.xyz << EOF
3

OW     2.30    6.28    1.13
HW     1.37    6.26    1.50 
HW     2.31    5.89    0.21
EOF

../yason.py water.yml | ../faunus --state state
if [[ $? == 0 ]]; then
    python ../yason.py < out.json > out.yml
    rm -f out.json
fi
