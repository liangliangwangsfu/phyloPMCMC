#!/usr/bin/env nextflow

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','')
// stuff to do

process summarizePipeline {

  cache false
  
  output:
      file 'pipeline-info.txt'
      
  publishDir deliverableDir, mode: 'copy', overwrite: true
  
  """
  echo 'scriptName: $workflow.scriptName' >> pipeline-info.txt
  echo 'start: $workflow.start' >> pipeline-info.txt
  echo 'runName: $workflow.runName' >> pipeline-info.txt
  echo 'nextflow.version: $workflow.nextflow.version' >> pipeline-info.txt
  """

}


process analysisCode {
cache true
  input:
    val gitRepoName from 'rejfreeAnalysis'
    val gitUser from 'alexandrebouchard'
    val codeRevision from 'deb3d12a7bfb361a694e21a357cf39ab4c9dbd12'
    val snapshotPath from '/home/shijia57/projects/def-liang-ab/shijia57/rejfreeAnalysis'
  
  output:
    file 'code' into analysisCode

  script:
    template 'buildSnapshot.sh'
}


process buildCode {
cache true
  input:
    val gitRepoName from 'phyloPMCMC'
    val gitUser from 'liangliangwangsfu'
    val codeRevision from '4ec4e1a96da1edd2648ba1df628d815b3cf67b62'
    val snapshotPath from '/home/shijia57/projects/def-liang-ab/shijia57/phyloPMCMC'
  
  output:
    file 'code' into code

  script:
    template 'buildSnapshot.sh' 
}



process run {

module 'java/1.8.0_121'  
echo true
  time '72h'
  memory '80 GB'
  cpus 2
  executor 'slurm'
  clusterOptions  '--account=def-liang-ab'

  input:
    file code
    
    each mainRand from  283..285
    each rand from  13199
    each particle from 5000
    each iterScalings from 5000
    echo true
        
  output:
    file '.' into execFolder
    
  """
  pwd  
  mkdir state
  mkdir state/execs
export LD_LIBRARY_PATH=$NIXUSER_PROFILE/lib:$LD_LIBRARY_PATH  
java   -cp 'code/lib/*' -Xmx80G   phyloPMCMC.PGSExperiments -useDataGenerator false    -nThousandIters 1 \
     -dataFile /home/shijia57/projects/def-liang-ab/shijia57/phyloPMCMC/Realdata/fish.msf  \
     -refTree /home/shijia57/projects/def-liang-ab/shijia57/phyloPMCMC/Realdata/consensus_fish.newick \
     -nTax 20  \
     -treeRate 10\
     -mainRand $mainRand  \
     -len  1000 \
     -sequenceType DNA \
     -generateDNAdata true \
     -nThreads 2 \
     -sampleTrans2tranv true\
     -saveTreesFromPMCMC false \
     -iterScalings  $iterScalings\
     -methods  PGS4K2P  \
     -useNonclock false \
     -useSlightNonclock false \
     -nParticlesEachStep $particle  \
     -nReplica 1 \
     -repPerDataPt  10 \
     -neighborPath  '/home/shijia57/bin/phylip-3.69/exe//neighbor' \
     -gen.rand $rand
     unset  LD_LIBRARY_PATH    
  """
}

process aggregate {

  input:
    file analysisCode
    file 'exec_*' from execFolder.toList()
    
  output:
    file aggregated    

  
  """
  ./code/bin/csv-aggregate \
    --experimentConfigs.managedExecutionFolder false \
    --experimentConfigs.saveStandardStreams false \
    --experimentConfigs.recordExecutionInfo false \
    --argumentFileName state/execs/0.exec/options.map \
    --argumentsKeys    gen.rand  PGSExperiments.iterScalings PGSExperiments.methods PGSExperiments.mainRand PGSExperiments.nParticlesEachStep\
    --dataPathInEachExecFolder state/execs/0.exec/results/logZout.csv \
    --outputFolderName aggregated
  """

}



process createPlot {

  echo true
  module 'r/3.4.0'
  input:
    file aggregated
    env SPARK_HOME from "${System.getProperty('user.home')}/bin/spark-2.1.0-bin-hadoop2.7"
        
   output:
 
    file 'logZ1.eps'
    file 'treeDistance.eps'


    
  publishDir deliverableDir, mode: 'copy', overwrite: true
  
  afterScript 'rm -r metastore_db; rm derby.log'
    
  """
  #!/usr/bin/env Rscript
  require("ggplot2")

  library(SparkR, lib.loc = c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib")))
  sparkR.session(master = "local[*]", sparkConfig = list(spark.driver.memory = "4g")) 
   
  data <- read.df("$aggregated", "csv", header="true", inferSchema="true")
  head(data)
  data <- collect(data)
  nR <- nrow(data)
  nrow <- which(data[,2] == 'MCMC')
  data2 <- data[-nrow,]


  gname = c("logZ1.eps",sep="")  
  postscript(gname,width=7,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
  par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

  p <- ggplot(data2, aes(Method, logZ))
  p + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = logZ))+ theme(legend.title=element_blank())+theme_bw()
  dev.off()
  
  require("Sleuth2")
  require("ggpubr")
  results <- read.df("$aggregated2", "csv", header="true", inferSchema="true")
  results <- collect(results)
  head(results)
  head(results[-which(results[,c("Metric")] == 'BestSampledLogLL'),])
  results2 <- data.frame(results[-which(results[,c("Metric")] == 'BestSampledLogLL'),])
  results2[,2] <- paste(results2[,1], results2[,2])
  head(results2)
  
  p <- ggplot(results2[which(results2[,7] == 'ConsensusLogLL'),], aes(Adaptive, Value))

  p2 <- ggplot(results2[which(results2[,7] == 'PartitionMetric'),], aes(Adaptive, Value))

  p3 <- ggplot(results2[which(results2[,7] == 'RobinsonFouldsMetric'),], aes(Adaptive, Value))

  p4 <- ggplot(results2[which(results2[,7] == 'KuhnerFelsenstein'),], aes(Adaptive, Value))

  gname = c("treeDistance.eps",sep="")  
  postscript(gname,width=8,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
  #par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

  ggarrange(p + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1)+rremove("x.text")+rremove("ylab")+rremove("legend")+ geom_boxplot(aes(color = Adaptive))+ xlab('ConsensusLogLL')
          ,p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1)+rremove("x.text")+rremove("ylab") + geom_boxplot(aes(color = Adaptive))+ xlab('PartitionMetric')
          ,p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1)+rremove("x.text")+rremove("ylab") + geom_boxplot(aes(color = Adaptive))+ xlab('RobinsonFouldsMetric')
          ,p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1)+rremove("x.text")+rremove("ylab") + geom_boxplot(aes(color = Adaptive))+ xlab('KuhnerFelsenstein')
          ,ncol = 4, nrow = 1, common.legend = TRUE)

  dev.off()
  getwd()

  """
}


