Summary
-------
phyloPMCMC is a novel  combinatorial  sequential Monte Carlo (CSMC) method by  proposing a more efficient proposal  distribution.
We combine the CSMC and MCMC in the framework of the particle Gibbs (PG) sampler to jointly estimate the phylogenetic trees and evolutionary parameters. The new algorithm can be easily parallelized by allocating samples over different computing cores.
We validate that the developed CSMC can sample trees more efficiently in various particle Gibbs samplers  via numerical experiments.




Installation
------------


There are several options available to install the package:

### Integrate to a gradle script

Simply add the following lines:

```groovy
repositories {
 mavenCentral()
 jcenter()
 maven {
    url "http://www.stat.sfu.ca/~lwa68/maven/"
  }
}

dependencies {
  compile group: 'ca.sfu.stat', name: 'phyloPMCMC', version: '1.0.0'
}
```

### Compile using the provided gradle script

- Check out the source ``git clone git@github.com:liangliangwangsfu/phyloPMCMC.git``
- Compile using ``./gradlew installDist``
- Add the jars in ``build/install/phyloPMCMC/lib/`` into your classpath

### Use in eclipse

- Check out the source ``git clone git@github.com:liangliangwangsfu/phyloPMCMC.git``
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root
  - Deselect ``Copy projects into workspace`` to avoid having duplicates


Usage
-----

### Quick start
PGSExperiments -useDataGenerator true
  -treeRate 10 -nThousandIters 1 -nTax 10 -len  100 -sequenceType DNA -generateDNAdata true   -useNonclock false -useSlightNonclock false -nThousandIters 1 -iterScalings  300  -methods   PGS4K2PBF -nCSMC 2  -nUCSMC 2 -resamplingStrategy ESS  -essRatioThreshold 0.5 -sampleTrans2tranv true -nReplica 1 -mainRand  289 -gen.rand  16   -saveTreesFromPMCMC true -nThreads 2  -nParticlesEachStep 300 -repPerDataPt  1 -mrBayesPath /Users/oudomame/Dropbox/phyloSoftware/mrbayes-3.2.6/src//mb  -neighborPath /Users/oudomame/Dropbox/phyloSoftware/phylip-3.69/exe//neighbor

For settings in our experimental results, please refers to folder ``setups''.

