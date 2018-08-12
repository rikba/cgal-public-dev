# Generate the benchmark table in the user manual

First, compile the file `benchmark_region_growing_with_points_2.cpp` using the provided cmake file:

```bash
$ cd /path/to/build/directory
$ cmake -DCGAL_DIR=/path/to/cgal/release/build -DCMAKE_BUILD_TYPE=Release /path/to/benchmark/Generalized_region_growing
$ make benchmark_region_growing_with_points_2 # Compile
```

The program uses data file `data/inputbig_2.xyz`, run the benchmark program as follows:

```bash
$ ./benchmark_region_growing_with_points_2 /path/to/data/inputbig_2.xyz
```

Result:

```
Test #1:
  radius = 1;
  min_size = 5;
  epsilon = 4.5;
  normal_threshold = 0.7;
  -----
  Time elapsed: 0.350881
  Number of regions detected: 796
  Number of points assigned: 4491
Test #2:
  radius = 3;
  min_size = 5;
  epsilon = 4.5;
  normal_threshold = 0.7;
  -----
  Time elapsed: 0.306144
  Number of regions detected: 3054
  Number of points assigned: 63154
Test #3:
  radius = 6;
  min_size = 5;
  epsilon = 4.5;
  normal_threshold = 0.7;
  -----
  Time elapsed: 0.435842
  Number of regions detected: 2483
  Number of points assigned: 64977
Test #4:
  radius = 9;
  min_size = 5;
  epsilon = 4.5;
  normal_threshold = 0.7;
  -----
  Time elapsed: 0.650403
  Number of regions detected: 2282
  Number of points assigned: 65353
```
