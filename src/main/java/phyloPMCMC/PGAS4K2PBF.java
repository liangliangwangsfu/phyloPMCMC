package phyloPMCMC;

import static nuts.util.CollUtils.list;
import java.io.File;
import java.util.List;
import java.util.Random;
import nuts.io.IO;
import nuts.math.Sampling;
import pty.RootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.mcmc.UnrootedTreeState;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleKernel;
import pty.smc.models.CTMC;
import smc.BackForwardKernel;
import smc.PartialCoalescentState4BackForwardKernel;
import smc.PGASParticleFilter;
import smc.PGASParticleFilter.ResamplingStrategy;
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
	private String nameOfAllTrees="allTrees-PGAS4K2PBF.trees";
	private boolean saveTreesFromPMCMC=false;
	private int sampleTreeEveryNIter=1; 
	private boolean processTree=false;
	private boolean isGS4Clock=true;
	private double trans2tranv=2;
	public double a=1.25;
	private PartialCoalescentState4BackForwardKernel sampled=null;
	private boolean sampleTrans2tranv=true;
	private CTMC ctmc =null;
	private PartialCoalescentState4BackForwardKernel init=null;

	public boolean isSampleTrans2tranv() {
		return sampleTrans2tranv;
	}

	public void setSampleTrans2tranv(boolean sampleTrans2tranv) {
		this.sampleTrans2tranv = sampleTrans2tranv;
	}

	public PGAS4K2PBF(Dataset dataset0,ParticleFilterOptions options,TreeDistancesProcessor tdp,
			boolean useTopologyProcessor,TreeTopologyProcessor trTopo,
			RootedTree initrt, boolean processTree,boolean isGS4Clock,int sampleTreeEveryNIter,PartialCoalescentState4BackForwardKernel init, CTMC ctmc)
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
		this.init=init;
		this.ctmc=ctmc;

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

	private double MHTrans2tranv(double currentTrans2tranv, Random rand) {
		Trans2tranvProposal kappaProposal=new Trans2tranvProposal(a,rand);
		Pair<Double,Double> proposed=kappaProposal.propose(currentTrans2tranv);		
		double proposedTrans2tranv=proposed.getFirst();
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);			
		UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);
		double logratio = ncs.logLikelihood() - previousLogLLEstimate+proposed.getSecond();
		double acceptPr = Math.min(1, Math.exp(logratio));
		final boolean accept = Sampling.sampleBern(acceptPr, rand);
		if (accept) {
			trans2tranv = proposedTrans2tranv;
			//	previousLogLLEstimate = ncs.logLikelihood();			
		}
		return acceptPr;
	}

	public static class Trans2tranvProposal{
		private final double a; // higher will have lower accept rate		
		private Random rand;
		public Trans2tranvProposal(double a, Random rand) {
			if (a <= 1)
				throw new RuntimeException();
			this.a = a;
			this.rand=rand;			
		}

		public  Pair<Double,Double> propose(double currentTrans2tranv) {
			double lambda=2*Math.log(a);
			double rvUnif = Sampling.nextDouble(rand, 0, 1);
			double m  = Math.exp(lambda*(rvUnif-0.5));		
			final double newTrans2tranv = m * currentTrans2tranv;
			return Pair.makePair(newTrans2tranv, Math.log(m));						
		}
	}

	public void next(Random rand)
	{
		iter++;
		RootedTree previousSample = currentSample;		
		if(sampleTrans2tranv) MHTrans2tranv(trans2tranv,  rand);
		// sample from PF
		StoreProcessor<PartialCoalescentState4BackForwardKernel> pro = new StoreProcessor<PartialCoalescentState4BackForwardKernel>();		 
		if((iter % sampleTreeEveryNIter) == 0)
		{
			if(sampleTrans2tranv) {
			ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), trans2tranv);
            PartialCoalescentState init0 = PartialCoalescentState.initFastState(dataset, ctmc, true);			
            init = new PartialCoalescentState4BackForwardKernel(init0, null, null, null,  0, new int[] {-1,-1});			 
			}
		
			
			ParticleKernel<PartialCoalescentState4BackForwardKernel> kernel = new BackForwardKernel(
					init);
			PGASParticleFilter<PartialCoalescentState4BackForwardKernel> pf = new PGASParticleFilter<PartialCoalescentState4BackForwardKernel>();
			pf.rand= rand;
			pf.nThreads = options.nThreads;
			pf.resampleLastRound = false;
			pf.N=options.nParticles;
			//pf.resamplingStrategy=ResamplingStrategy.ESS;
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
			currentSample=sampled.getCurrentState().getFullCoalescentState();
			//			previousLogLLEstimate=sampled.logLikelihood();
			UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);
			previousLogLLEstimate=ncs.logLikelihood();  //TODO:update the logLikelihood calculation in PartialCoalescentState4BackForwardKernel so that it is equal to this value.
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
				"trans2tranv", trans2tranv,
				//        "maskSparsity", currentSparsity,
				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
				"LogLikelihood", previousLogLLEstimate);
	}

}



