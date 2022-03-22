## Granger Components Analysis Homepage

GCA is an unsupervised learning algorithm that automatically identifies Granger causal connections in your data!

Here you can find the latest implementations, beginning with the original Matlab code.

```markdown
[W,V] = runGca(X)
```

### load the provided simulated data
```markdown
data=load('rez_simulated_var_2pair_15-Mar-2022.mat');
```

### run the algorithm on a single realization of a VAR(3) process
```markdown
[W,V,G,stats]=runGca(data.figStats.X,3,2);
```

### display results
```markdown
[W,V,G,stats]=runGca(data.figStats.X,3,2);
```

[comment]: [<img src="https://github.com/dmochow/gca/blob/main/arxiv-logo-1.png">](https://arxiv.org/abs/2203.10679)

[comment]: You can use the [editor on GitHub](https://github.com/dmochow/gca/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.




