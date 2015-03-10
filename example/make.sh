#!/bin/bash
( cd convex && cmake . -DCMAKE_BUILD_TYPE=Release && make )
( cd analysis && cmake . -DCMAKE_BUILD_TYPE=Release && make )
( cd ball_knob && EXACT=1 cmake . -DCMAKE_BUILD_TYPE=Release && make )
( cd find_path && EXACT=1 cmake . -DCMAKE_BUILD_TYPE=Release && make )
( cd gearbox && cmake . -DCMAKE_BUILD_TYPE=Release && make )
( cd simple_raster && EXACT=1 cmake . -DCMAKE_BUILD_TYPE=Release && make )
