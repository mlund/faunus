#!/bin/bash

# create water.xyz
cat > water.xyz << EOF
3

OW     2.30    6.28    1.13
HW     1.37    6.26    1.50 
HW     2.31    5.89    0.21
EOF

python yason.py < water.yml > water.json # yml -> json
rc=$?
if [[ $rc != 0 ]]; then
    echo "error parsing yaml file"
    exit $rc
fi

if [ -f example ]; then
    ./example < water.json 
    if [[ $? == 0 ]]; then
        python yason.py < out.json > out.yml
        rm out.json
    fi
fi
