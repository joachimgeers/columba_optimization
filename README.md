# columba_optimization
Sequential optimization framework for Columba (https://github.com/biointec/columba). Uploaded for my master's thesis.

Simply add the src files to Columba's src files, and replace the CMakeFile. Next you can follow the compilation steps as given in Columba's repository (note: the optimization framework only works for custom search schemes).

The optimization framework can then be run as follows:

```columba_optimize_weight -ss custom "path_to_custom_scheme" "[basefile] [readfile]```

Options:

```
--seed                  The seeding value for sampling.
--threshold             Value for the sampling threshold
--reheat-factor         The restart factor for the framework
--reheat-iterations     The number of iterations before reheating
--accepts               The number of accepts before changing the SCHC threshold
--max-weight            The maximum weight parameter for the dynamic partitioning while training
--sample-size           The size of the training sample
--duration              The length of training
-e                      Distance
-m                      Distance metric
-p                      Partitioning technique
```
