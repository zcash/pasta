#!/bin/sh
# pip install bintrees Pillow
# apt-get install ffmpeg ffcvt

set -e
./checksumsets.py --animate
if [ -f animation-p.gif ]; then
    ffcvt -f animation-p.gif
    mv -f animation-p_.webm animation-p.webm
    rm -f animation-p.gif
fi
if [ -f animation-q.gif ]; then
    ffcvt -f animation-q.gif
    mv -f animation-q_.webm animation-q.webm
    rm -f animation-q.gif
fi
