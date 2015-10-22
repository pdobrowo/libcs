#!/bin/bash
( cd convex && cmake . -DCMAKE_BUILD_TYPE=Debug && make )
( cd analysis && cmake . -DCMAKE_BUILD_TYPE=Debug && make )
( cd ball_knob && EXACT=1 cmake . -DCMAKE_BUILD_TYPE=Debug && make )
( cd find_path && EXACT=1 cmake . -DCMAKE_BUILD_TYPE=Debug && make )
( cd gearbox && cmake . -DCMAKE_BUILD_TYPE=Debug && make )
( cd simple_raster && EXACT=1 cmake . -DCMAKE_BUILD_TYPE=Debug && make )
