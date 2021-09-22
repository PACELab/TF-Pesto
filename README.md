# TF-Pesto
A Distributed Version of TensorFlow r1.14 to support User-specified Operation Partitioning and Scheduling 

## Docker Image:
https://hub.docker.com/repository/docker/xilenteyex/tf_pesto

## Requirements:
Docker with Nvidia-Runtime

## How to Run:

```
$ docker pull xilenteyex/tf_pesto:latest
$ docker run -it --runtime=nvidia xilenteyex/tf_pesto:latest bash
```

## Custom Device Placement and Scheduling:
```
$ cd /root/
$ mkdir pestoPlacement
$ cd pestoPlacement
```


```
$ vi part.json
```
Populate this file with the partitions in following json format:

``` json
[
  {"op-Name":"requested-Device"},
  {"op-Name":"requested-Device"},
  ...
  {"op-Name":"requested-Device"}
]
```

```
$ vi cdep.json
```

Populate this file with the Scheduling constraints in following json format:

``` json
[
  {"op-Name":["<predecessor-OpName-list>"]},
  {"op-Name":["<predecessor-OpName-list>"]},
  ...
  {"op-Name":["<predecessor-OpName-list>"]}
]
```

## Build TensorFlow
```
$ cd /tensorflow_src/
$ bazel build --config=opt --config=cuda //tensorflow/tools/pip_package:build_pip_package
$ ./bazel-bin/tensorflow/tools/pip_package/build_pip_package /tmp/
$ pip uninstall -y tensorflow
$ pip install /tmp/tensorflow-1.14.0-cp27-cp27mu-linux_x86_64.whl
```

#### Finally, we can run the TensorFlow program as usual and it will follow the placements specified in above format


## Simulator

### Requirements:
Matlab with Bioinformatics Toolbox and Cplex

### How to Run:
Simply execute main.m
```
Configurable parameters:
timelmt: Time limit for the solver;
para_config: [maximum size of the supernode dag; maximum size of all supernodes; maximum length of all supernodes];
speed_mat : 2*2 matrix storing parameters from linear interpolation of CPU-GPU speed and GPU-GPU speed.
```
