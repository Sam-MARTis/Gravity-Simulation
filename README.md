# N body gravity simulator


The main point of this repo is to explore the barnes hut algorithm in more detail and integrate it with the GPU to simulate a large number of particles.

To run, go to the correct folder and run the `grav` executable
### Barnes Hut Algorithm
The conventional naive algorithm for n-body sim involves calculating the effect of each particle on the particle in question and then repeat that for every particle. For N particles, that is of the time complexity of O(N^2), which isn't great. The BH algorithm is an O(NlogN) algorithm, which is much better. It's based on the fast multipole method. 
In essence, if a group of particles are huddled close together OR the group is quite far away, the group can be treated as a single particle with the mass equal to the group's mass at the COM location. 

Parallelizing this is incredibly hard cause of the irregular nature of the Quadtree/Octtree used to partition space in the implementation of BH.

Based on the following paper:
https://iss.oden.utexas.edu/Publications/Papers/burtscher11.pdf
An Efficient CUDA Implementation of the Tree-Based Barnes Hut n-Body Algorithm

Two main milestones:
- CPU implementation
- GPU implementation (TODO)

So heres what I'm thinking of doing in the below order:
- CPU naive algorithm
- CPU barnes hut
- CPU naive algroithm multithreading
- CPU barnes hut multithreading
- GPU naive algorithm
- GPU barnes hut

CPU naive O(N^2) algorithm
![01](01.png)

CPU barnes hut O(NlogN) algorithm + Collision detection
![02](02.png)
CPU naive algorithm multithreading
![03](03.png)
CPU barnes hut multithreading
![04](04.png)
05- Todo