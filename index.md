## Granger Components Analysis Homepage

GCA is an unsupervised learning algorithm that automatically identifies Granger causal connections in your data!

Here you can find the latest implementations, beginning with the original Matlab code.

### Load the provided simulated data
```markdown
data=load('rez_simulated_var_2pair_15-Mar-2022.mat');
```

### Run the algorithm on a single realization of a VAR(3) process
```markdown
[W,V,G,stats]=runGca(data.figStats.X,3,2);
```

### Display results to show convergence for first two component pairs
```markdown
figure; plot(stats.fvals_t); ylabel('-G'); legend('Pair 1', 'Pair 2');
```

![demo result](https://github.com/dmochow/gca/blob/main/demo_result.png)

![<img src="https://github.com/dmochow/gca/blob/main/arxiv-logo-1.png">](https://arxiv.org/abs/2203.10679)





