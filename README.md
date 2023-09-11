# 2nd-DTM: Direct Transcription method for Second-order Systems
2nd-DTM contains numerical methods to solve continuous-time trajectory optimization for second-order systems. I developed it while reading a very helpful [tutorial](https://www.matthewpeterkelly.com/tutorials/trajectoryOptimization/index.html) on trajectory optimization by Kelly. The first-order methods in 2nd-DMT are mainly based on [Kelly's implementation](https://github.com/MatthewPeterKelly/OptimTraj). In addition, I developed some methods for systems with second-order dynamics. Check it by running the `examples` in this repo.

The code is powered by CasADi.


## Usage
1. download casadi from [here](https://web.casadi.org/get/)
2. Modify addpath on casadi in setup.m
3. Run setup.m
4. Run demo.m

Check folder `examples` for more examples.

## Learn more about trajectory optimimization
1. Very helpful [tutorial](https://www.matthewpeterkelly.com/tutorials/trajectoryOptimization/cartPoleCollocation.svg), [video](https://www.youtube.com/watch?v=wlkRYMVUZTs&feature=youtu.be), [code](https://github.com/MatthewPeterKelly/OptimTraj), and [conresponding paper](https://epubs.siam.org/doi/10.1137/16M1062569) by Kelly.
2. Our [2nd-Collocation]() and [2nd-Shooting]() papers (Will be available later).
   
