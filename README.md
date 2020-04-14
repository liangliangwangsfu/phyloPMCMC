Summary
-------

<!-- [![Build Status](https://travis-ci.org/alexandrebouchard/phyloPMCMC.png?branch=master)](https://travis-ci.org/alexandrebouchard/phyloPMCMC) -->

phyloPMCMC is a novel  combinatorial  sequential Monte Carlo (CSMC) method by  proposing a more efficient proposal  distribution.
We combine the CSMC and MCMC in the framework of the particle Gibbs (PG) sampler to jointly estimate the phylogenetic trees and evolutionary parameters. The new algorithm can be easily parallelized by allocating samples over different computing cores.
We validate that the developed CSMC can sample trees more efficiently in various particle Gibbs samplers  via numerical experiments.




Installation
------------


There are several options available to install the package:

### Integrate to a gradle script

Simply add the following lines (replacing 1.0.0 by the current version (see git tags)):

```groovy
repositories {
 mavenCentral()
 jcenter()
 maven {
    url "http://www.stat.ubc.ca/~bouchard/maven/"
  }
}

dependencies {
  compile group: 'ca.ubc.stat', name: 'phyloPMCMC', version: '1.0.0'
}
```

### Compile using the provided gradle script

- Check out the source ``git clone git@github.com:alexandrebouchard/phyloPMCMC.git``
- Compile using ``./gradlew installDist``
- Add the jars in ``build/install/phyloPMCMC/lib/`` into your classpath

### Use in eclipse

- Check out the source ``git clone git@github.com:alexandrebouchard/phyloPMCMC.git``
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root
  - Deselect ``Copy projects into workspace`` to avoid having duplicates


Usage
-----

### Quick start

...
