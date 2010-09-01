#!/bin/sh -e

set -e

(cd qt && sh build.sh)
(cd dashboard && sh build.sh)
(cd swing && sh build.sh)

echo "^C to stop or Enter to upload"
read
sh upload.sh
