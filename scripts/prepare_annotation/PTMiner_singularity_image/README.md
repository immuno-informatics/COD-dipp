# Step 1: build the image on a local PC with sudo access
```
sudo singularity build --sandbox ubuntu_wine.simg ubuntu20.04_wine.singularity
```

# Step 2: launch the image to initialize it on the local PC
```
singularity exec -w ubuntu_wine.simg /init
```

# Step 3: convert the image to read only on the local PC
```
sudo singularity build PTMiner.simg ubuntu_wine.simg/
```

# Step 4: send the read only image from the local PC to the HPC


# example run 
```
# wine can find "/DATA" using "z:\DATA" path
# test.param: PTMiner parameter file
singularity exec -B $PWD/test_data:/data ./PTMiner.simg /PTMiner /data/test.param
```

# launch on test data example
```
singularity exec -B $PWD/test_data:/data ./PTMiner.simg /PTMiner /data/test.param
```

