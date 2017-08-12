package phyloPMCMC;

import static nuts.util.CollUtils.list;
import java.io.File;
import java.util.List;
import java.util.Random;
import nuts.io.IO;
import pty.RootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleKernel;
import pty.smc.models.CTMC;
import smc.BackForwardKernel;
import smc.PartialCoalescentState4BackForwardKernel;
import smc.PGASParticleFilter;
import smc.PGASParticleFilter.StoreProcessor;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import fig.basic.Pair;
import fig.exec.Execution;
import gep.util.OutputManager;

public class PGAS4K2PBF {
	private final Dataset dataset;
	ParticleFilterOptions options=null;
	private final TreeDistancesProcessor tdp;
	private boolean useTopologyProcessor=false;	
	private final TreeTopologyProcessor trTopo;
	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
	private RootedTree currentSample = null;
	public static OutputManager outMan = new OutputManager();		
	private int iter=0; 
	private int treeCount=0;	
	File output=new File(Execution.getFile("results")); 
	private String nameOfAllTrees="allTrees.trees";
	private boolean saveTreesFromPMCMC=false;
	private int sampleTreeEveryNIter=1; 
	private boolean processTree=false;
	private boolean isGS4Clock=true;
	private double trans2tranv;
	public double a=2;
	private PartialCoalescentState4BackForwardKernel sampled=null;


	public PGAS4K2PBF(Dataset dataset0,ParticleFilterOptions options,TreeDistancesProcessor tdp,
			boolean useTopologyProcessor,TreeTopologyProcessor trTopo,
			RootedTree initrt, boolean processTree,boolean isGS4Clock,int sampleTreeEveryNIter)
	{
		this.dataset=dataset0;
		this.options = options;
		this.tdp = tdp;
		this.useTopologyProcessor=useTopologyProcessor;
		if(useTopologyProcessor)
			this.trTopo=trTopo; 
		else this.trTopo=null; 
		this.currentSample=initrt;   
		this.processTree=processTree;
		this.isGS4Clock=isGS4Clock;
		this.sampleTreeEveryNIter=sampleTreeEveryNIter;
	}

	public PGAS4K2PBF(Dataset dataset0, ParticleFilterOptions options,  TreeDistancesProcessor tdp, RootedTree initrt, TreeTopologyProcessor trTopo)
	{
		this.dataset=dataset0;
		this.options = options;
		this.tdp = tdp;
		this.currentSample=initrt;
		this.trTopo=trTopo;
	}	

	public void setSaveTreesFromPMCMC(boolean saveTreesFromPMCMC)
	{
		this.saveTreesFromPMCMC=saveTreesFromPMCMC; 
	}

	public void setNameOfAllTrees(String nameOfAllTrees)
	{
		this.nameOfAllTrees=nameOfAllTrees; 
	}

	public String getNameOfAllTrees()
	{
		return this.nameOfAllTrees;
	}

	public boolean getSaveTreesFromPMCMC()
	{
		return this.saveTreesFromPMCMC;
	}

	public void setProcessTree(boolean processTree)
	{
		this.processTree=processTree; 
	}

	public RootedTree getRootedTree()
	{
		return currentSample;
	}

	public void next(Random rand)
	{
		iter++;
		RootedTree previousSample = currentSample;		
		// proposals:
		// alpha: multiplier				    
		//	double scale=Sampling.nextDouble(rand, 1.0/a, a);
		//		System.out.println(scale);
		//	double proposedTrans2tranv = scale*trans2tranv; 
		double proposedTrans2tranv = 2.0;
		// sample from PF
		StoreProcessor<PartialCoalescentState4BackForwardKernel> pro = new StoreProcessor<PartialCoalescentState4BackForwardKernel>();		 
		if((iter % sampleTreeEveryNIter) == 0)
		{
			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);
			PartialCoalescentState init0 = PartialCoalescentState
					.initFastState(dataset, ctmc, true);
			PartialCoalescentState4BackForwardKernel init = new PartialCoalescentState4BackForwardKernel(
					init0, null, 0, new int[] {-1,-1});

			ParticleKernel<PartialCoalescentState4BackForwardKernel> kernel = new BackForwardKernel(
					init);

			PGASParticleFilter<PartialCoalescentState4BackForwardKernel> pf = new PGASParticleFilter<PartialCoalescentState4BackForwardKernel>();
			pf.rand= rand;
			pf.nThreads = options.nThreads;
			pf.resampleLastRound = false;
			pf.N=options.nParticles;
			if(sampled!=null)
			{
			//	System.out.println("Find the conditioned path!");
			List<Pair<PartialCoalescentState4BackForwardKernel, Double>> restorePCS = PartialCoalescentState4BackForwardKernel.restoreSequence(sampled);			
			List<PartialCoalescentState4BackForwardKernel> path = list();
			double[] weights=new double[restorePCS.size()];
			for(int i=0;i<restorePCS.size();i++){
				path.add(restorePCS.get(i).getFirst());			 
				weights[i]=restorePCS.get(i).getSecond();
			}
			// set the conditioning and its weights
			pf.setConditional(path, weights);
			}
			// do the sampling			
			pf.sample(kernel,  pro);
			sampled = pro.sample(rand);			
			currentSample=sampled.getFullCoalescentState();
			// update tdp
			if(processTree)tdp.process(currentSample);
			if(useTopologyProcessor) trTopo.process(currentSample);
			treeCount++;
			if(saveTreesFromPMCMC)
			{
				String stringOfTree=RootedTree.Util.toNewick(currentSample);
				String cmdStr="echo -n 'TREE gsTree_" + treeCount+"=' >>"+nameOfAllTrees;
				IO.call("bash -s",cmdStr,output);
				cmdStr="echo '" +stringOfTree + "' | sed 's/internal_[0-9]*_[0-9]*//g'"+ " | sed 's/leaf_//g'"+ " >> "+nameOfAllTrees;
				IO.call("bash -s",cmdStr,output);
			}
		}
		// log some stats
		final int tSize = currentSample.topology().nLeaves();
		outMan.write("PGAS4K2PBF",
				"Iter", iter,
				"treeSize", tSize,
				//        "maskSparsity", currentSparsity,
				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
				"LogLikelihood", previousLogLLEstimate);
	}

}



