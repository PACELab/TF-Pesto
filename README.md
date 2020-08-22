# TF-Pesto
A Distributed Version of TensorFlow r1.14 to support User-specified Operation Partitioning and Scheduling 


# Requirements:
Docker with Nvidia-Runtime

# How to Run:

```
$ docker pull xilenteyex/tf_pesto:latest
$ docker run -it --runtime=nvidia xilenteyex/tf_pesto:latest bash
```

# Custom Device Placement:
```
$ cd /root/
$ mkdir pestoPlacement
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
Populate this file with the scheduling constraints in json format:
[
  {"op-Name":["<predecessor-OpName-list>"]},
  {"op-Name":["<predecessor-OpName-list>"]},
  ...
  {"op-Name":["<predecessor-OpName-list>"]}
]
```
