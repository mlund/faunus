#!/usr/bin/env python

# This will download select Commodore 64 SID files from the
# online High Voltage SID Collection (HVSC) and generate a JSON playlist
# to be read by the faunus exe.

import sys
if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

import os, json, urllib.request
try:
    import ruamel.yaml as yaml
except ImportError:
    import yaml

if len(sys.argv)!=2:
    print("usage: {} music.yml".format(sys.argv[0]))
    sys.exit(1)

filename = sys.argv[1]

with open(filename) as f:
    sidlist = yaml.safe_load(f)
    server = sidlist['server'] + '/'
    dstdir = 'sids/'
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    for i in sidlist['songs']:
        url = server + i['file']
        basename = os.path.basename(url)
        dstfile = dstdir + basename
        i['file'] = dstfile
        if not os.path.isfile(dstfile):
            print('retrieving', dstfile)
            urllib.request.urlretrieve(url, dstfile)

    # save playlist to json
    out = json.dumps(sidlist, indent=2)
    with open(dstdir + 'music.json', 'w') as f:
        f.write(out)
