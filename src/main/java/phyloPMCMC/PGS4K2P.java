package phyloPMCMC;

import static nuts.util.CollUtils.list;
import static nuts.util.CollUtils.map;
import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import conifer.trees.StandardNonClockPriorDensity;
import nuts.io.IO;
import nuts.math.Sampling;
import nuts.util.Arbre;
import pty.RootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.mcmc.UnrootedTreeState;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.models.CTMC;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import fig.basic.Pair;
import fig.exec.Execution;
import fig.prob.Gamma;
import gep.util.OutputManager;
import goblin.Taxon;

public class PGS4K2P {
	private final Dataset dataset;
	ParticleFilterOptions options = null;
	private final TreeDistancesProcessor tdp;
	private boolean useTopologyProcessor = false;
	private final TreeTopologyProcessor trTopo;
	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
	private RootedTree currentSample = null;
	public static OutputManager outMan = new OutputManager();
	private int iter = 0;
	private int treeCount = 0;
	File output = new File(Execution.getFile("results"));
	private String nameOfAllTrees = "allTrees-PGS4K2P.trees";
	private boolean saveTreesFromPMCMC = false;
	private int sampleTreeEveryNIter = 100;
	private boolean processTree = false;
	private boolean isGS4Clock = true;
	private double trans2tranv=2;
	public double a = 1.25;
    private boolean sampleTrans2tranv=true;
	
	public PGS4K2P(Dataset dataset0, ParticleFilterOptions options, TreeDistancesProcessor tdp,
			boolean useTopologyProcessor, TreeTopologyProcessor trTopo, RootedTree initrt, boolean processTree,
			boolean isGS4Clock, int sampleTreeEveryNIter) {
		this.dataset = dataset0;
		this.options = options;
		this.tdp = tdp;
		this.useTopologyProcessor = useTopologyProcessor;
		if (useTopologyProcessor)
			this.trTopo = trTopo;
		else
			this.trTopo = null;
		this.currentSample = initrt;
		this.processTree = processTree;
		this.isGS4Clock = isGS4Clock;
		this.sampleTreeEveryNIter = sampleTreeEveryNIter;
	}

	public PGS4K2P(Dataset dataset0, ParticleFilterOptions options, TreeDistancesProcessor tdp, RootedTree initrt,
			TreeTopologyProcessor trTopo) {
		this.dataset = dataset0;
		this.options = options;
		this.tdp = tdp;
		this.currentSample = initrt;
		this.trTopo = trTopo;
	}

	public void setSaveTreesFromPMCMC(boolean saveTreesFromPMCMC) {
		this.saveTreesFromPMCMC = saveTreesFromPMCMC;
	}

	public void setNameOfAllTrees(String nameOfAllTrees) {
		this.nameOfAllTrees = nameOfAllTrees;
	}

	public String getNameOfAllTrees() {
		return this.nameOfAllTrees;
	}

	public boolean getSaveTreesFromPMCMC() {
		return this.saveTreesFromPMCMC;
	}

	public void setProcessTree(boolean processTree) {
		this.processTree = processTree;
	}

	public RootedTree getRootedTree() {
		return currentSample;
	}

	private double MHTrans2tranv(double currentTrans2tranv, Random rand) {
		Trans2tranvProposal kappaProposal=new Trans2tranvProposal(a,rand);
		Pair<Double,Double> proposed=kappaProposal.propose(currentTrans2tranv);	
		double proposedTrans2tranv=proposed.getFirst();		
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);	

		UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);	
		//bug in the next line...
		System.out.println("..............."+ncs.logLikelihood());
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
		System.out.println("..............."+sampleTrans2tranv);
		if(sampleTrans2tranv) MHTrans2tranv(trans2tranv,  rand);
		// sample from PF
		StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();		 
		if((iter % sampleTreeEveryNIter) == 0)
		{
			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), trans2tranv);  
			PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, ctmc, true);  // is clock	          								
			LazyParticleKernel kernel = new PriorPriorKernel(init);	
			ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();	
			pf.N=options.nParticles;
			pf.rand = rand;
			pf.nThreads = options.nThreads;
			List<Pair<PartialCoalescentState,Double>> restorePCS = restoreSequence(kernel, currentSample,isGS4Clock); 
			List<PartialCoalescentState> path=list();
			double[] weights=new double[restorePCS.size()];
			for(int i=0;i<restorePCS.size();i++){
				path.add(restorePCS.get(i).getFirst());			 
				weights[i]=restorePCS.get(i).getSecond();
			}			 
			// set the conditioning and its weights
			pf.setConditional(path, weights);
			// do the sampling			
			pf.sample(kernel, pro);
			PartialCoalescentState sampled = pro.sample(rand);
			currentSample=sampled.getFullCoalescentState();
			previousLogLLEstimate=sampled.logLikelihood();
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
		outMan.write("PGS4K2P",
				"Iter", iter,
				"treeSize", tSize,
				"trans2tranv", trans2tranv,
				//        "maskSparsity", currentSparsity,
				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
				"LogLikelihood", previousLogLLEstimate);
	}

	public static double height(Map<Taxon,Double> branchLengths, Arbre<Taxon> arbre)
	{
		double sum=0; 
		do{
			arbre=arbre.getChildren().get(0);
			sum+=branchLengths.get(arbre.getContents());         			   
		}while(!arbre.isLeaf());
		return sum;
	}

	public static Map sortByValue(Map map) {
		List list = new LinkedList(map.entrySet());
		Collections.sort(list, new Comparator() {
			public int compare(Object o1, Object o2) {
				return ((Comparable) ((Map.Entry) (o1)).getValue())
						.compareTo(((Map.Entry) (o2)).getValue());
			}
		});
		Map result = new LinkedHashMap();
		for (Iterator it = list.iterator(); it.hasNext();) {
			Map.Entry entry = (Map.Entry)it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	} 

	public static List<Pair<PartialCoalescentState,Double>> restoreSequence4NonClockTree(PartialCoalescentState current,RootedTree rt)
	{
		List<Pair<PartialCoalescentState,Double>> result=list();			 		
		List<Arbre<Taxon>> childrenList=rt.topology().nodes();
		class ArbreComparator implements Comparator<Arbre<Taxon>> {
			@Override
			public int compare(Arbre<Taxon> arbre1, Arbre<Taxon> arbre2) {		    	
				return   arbre1.getContents().toString().compareTo(arbre2.getContents().toString());      
			}
		}
		Collections.sort(childrenList, new ArbreComparator());
		Map<Taxon,Double> branchLengths=rt.branchLengths();
		for(int i=0;i<childrenList.size();i++)
		{
			if(!childrenList.get(i).isLeaf())
			{
				Arbre<Taxon> currentArbre=childrenList.get(i);
				Taxon first=currentArbre.getChildren().get(0).getContents(),second=currentArbre.getChildren().get(1).getContents();
				PartialCoalescentState coalesceResult = current.coalesce(current.indexOf(first), current.indexOf(second), 0, branchLengths.get(first), branchLengths.get(second),currentArbre.getContents());							
				double logWeight = coalesceResult.logLikelihoodRatio()-Math.log(coalesceResult.nNonTrivialRoots());  //coalesceResult.logLikelihood()-current.logLikelihood();								
				current=coalesceResult;
				result.add(Pair.makePair(coalesceResult, logWeight));
			}
		}		
		return result;
	}

	public static List<Pair<PartialCoalescentState,Double>> restoreSequence(ParticleKernel<PartialCoalescentState> kernel,RootedTree rt,boolean isClock)
	{
		if(!isClock)return restoreSequence4NonClockTree(kernel.getInitial(),rt);
		List<String> newNodeNames=list();
		List<Pair<PartialCoalescentState,Double>> result=list();
		PartialCoalescentState current= kernel.getInitial();
		List<Arbre<Taxon>> childrenList=rt.topology().nodes();
		Map<Taxon,Double> branchLengths=rt.branchLengths();
		Map<Arbre<Taxon>,Double> heightMap=map();
		for(int i=0;i<childrenList.size();i++)
		{
			if(!childrenList.get(i).isLeaf()) 
				heightMap.put(childrenList.get(i), height(branchLengths, childrenList.get(i))); 		
		}
		//   sort the height in an increasing order		
		Map<Arbre<Taxon>,Double> sortHeightMap=sortByValue(heightMap);
		Set<Arbre<Taxon>> arbreSet=sortHeightMap.keySet();		
		Iterator<Arbre<Taxon>> arbreIterator=arbreSet.iterator();
		while(arbreIterator.hasNext())
		{			
			Arbre<Taxon> currentArbre=arbreIterator.next();
			newNodeNames.add(currentArbre.toString());
		}
		arbreIterator=arbreSet.iterator();
		double previousHeight=0;
		while(arbreIterator.hasNext())
		{			
			Arbre<Taxon> currentArbre=arbreIterator.next();
			Taxon first=currentArbre.getChildren().get(0).getContents(),second=currentArbre.getChildren().get(1).getContents();
			double currentHeight=sortHeightMap.get(currentArbre);
			double currentDelta=currentHeight-previousHeight;
			previousHeight=currentHeight;
			PartialCoalescentState coalesceResult = current.coalesce(current.indexOf(first), current.indexOf(second), currentDelta, 0, 0,currentArbre.getContents());							
			double logWeight = coalesceResult.logLikelihoodRatio();// coalesceResult.logLikelihood()-current.logLikelihood();
			current=coalesceResult;
			result.add(Pair.makePair(coalesceResult, logWeight));			
		}			
		return result;
	}

	public boolean isSampleTrans2tranv() {
		return sampleTrans2tranv;
	}

	public void setSampleTrans2tranv(boolean sampleTrans2tranv) {
		this.sampleTrans2tranv = sampleTrans2tranv;
	}
}
