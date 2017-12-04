#!/bin/bash

# Generate two single atom molecules, A.xyz and B.xyz
cat > A.xyz << EOF
1

A      0.000   0.000  0.000
EOF
cat A.xyz | tr "A" "B" > B.xyz

python yason.py < example.yml > example.json # yml -> json
rc=$?
if [[ $rc != 0 ]]; then
    echo "error parsing yaml file"
    exit $rc
fi

if [ -f example ]; then
    ./example
    python yason.py < out.json
fi
