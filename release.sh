#!/bin/sh -e

set -e

(cd qt && sh build.sh)
(cd dashboard && sh build.sh)
(cd swing && sh build.sh)

echo "^C to stop or ^D to upload"
cat
sh upload.sh
