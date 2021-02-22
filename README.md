# Semi-Supervised Off Policy Reinforcement Learning

Code for Semi-Supervised Off Policy Reinforcement Learning. If you use our code please cite our [paper](https://arxiv.org/abs/2012.04809).

Repo is set up for Simulations shown in the manuscript and will work for any 2-time point, binary action setting. 

### Overview

Running `MVsims.R` will 

1) Generate data & get into suitable format, alternatively you can use your own data
2) Estimate a policy using supervised and semi-supervised Q-learning
3) Print Q-learning parameters and standard errors
3) Use doubly robust off policy evaluation for the estimated policies, implementation includes both supervised and semi-supervised
4) Print policy estimates and standard errors 

Results are saved in a list and saved in parallel directory `Results`

### Bibtex

```
@misc{sonabendw2021semisupervised,
      title={Semi-Supervised Off Policy Reinforcement Learning}, 
      author={Aaron Sonabend-W and Nilanjana Laha and Tianxi Cai and Rajarshi Mukherjee},
      year={2021},
      eprint={2012.04809},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```