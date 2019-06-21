package phyloPMCMC;

import static nuts.util.CollUtils.list;
import static nuts.util.CollUtils.map;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import nuts.io.IO;
import nuts.math.Sampling;
import nuts.util.Arbre;
import phyloPMCMC.PGS4K2PBF.Trans2tranvProposal;
import pty.RootedTree;
import pty.io.Dataset;
import pty.mcmc.UnrootedTreeState;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleFilter.ResamplingStrategy;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.models.CTMC;
import smc.BackForwardKernel;
import smc.PartialCoalescentState4BackForwardKernel;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import fig.basic.ListUtils;
import fig.basic.NumUtils;
import fig.basic.Pair;
import fig.exec.Execution;
import fig.prob.Dirichlet;
import gep.util.OutputManager;
import goblin.Taxon;

public class InteractingParticleGibbs4K2PBF {
	private final Dataset dataset;
	ParticleFilterOptions options = null;
	private final TreeDistancesProcessor tdp;
	private boolean useTopologyProcessor = false;
	private final TreeTopologyProcessor trTopo;
	private double[] previousLogLLEstimate; // = Double.NEGATIVE_INFINITY;
	private double LogLL;
	private List<RootedTree> currentSample = null;
	private double trans2tranv=2;
	public double a=1.25;
	private boolean sampleTrans2tranv=true;
	public static OutputManager outMan = new OutputManager();
	private int iter = 0;
	private int nCategories = 4;
	private int treeCount = 0;
	File output = new File(Execution.getFile("results"));
	private String nameOfAllTrees = "allTrees.trees";
	private boolean saveTreesFromPMCMC = false;
	private int sampleTreeEveryNIter = 100;
	private boolean processTree = false;
	private boolean isGS4Clock = true;
	private int nCSMC = 4; // the number of conditional SMC.
	private int nUCSMC = 4;// the number of unconditional SMC.
	private PartialCoalescentState4BackForwardKernel conditionedSampled=null;
	
	public boolean isSampleTrans2tranv() {
		return sampleTrans2tranv;
	}

	public void setSampleTrans2tranv(boolean sampleTrans2tranv) {
		this.sampleTrans2tranv = sampleTrans2tranv;
	}
	
	public InteractingParticleGibbs4K2PBF(Dataset dataset0,
			ParticleFilterOptions options, TreeDistancesProcessor tdp,
			boolean useTopologyProcessor, TreeTopologyProcessor trTopo,
			RootedTree initrt, 
			boolean processTree, boolean isGS4Clock, int sampleTreeEveryNIter,
			int nCSMC, int nUCSMC) {
		this.nCSMC = nCSMC;
		this.nUCSMC = nUCSMC;
		this.dataset = dataset0;
		this.options = options;
		this.tdp = tdp;
		this.useTopologyProcessor = useTopologyProcessor;
		if (useTopologyProcessor)
			this.trTopo = trTopo;
		else
			this.trTopo = null;
		this.currentSample = new ArrayList<RootedTree>(this.nCSMC);
		for (int i = 0; i < this.nCSMC; i++)
			this.currentSample.add(initrt);
		this.nCategories = nCategories;
		this.processTree = processTree;
		this.isGS4Clock = isGS4Clock;
		this.sampleTreeEveryNIter = sampleTreeEveryNIter;
		this.previousLogLLEstimate = new double[this.nCSMC];
		for (int i = 0; i < this.nCSMC; i++)
			previousLogLLEstimate[i] = Double.NEGATIVE_INFINITY;
	}
	
	public InteractingParticleGibbs4K2PBF(Dataset dataset0,
			ParticleFilterOptions options, TreeDistancesProcessor tdp,
			RootedTree initrt, TreeTopologyProcessor trTopo) { 
		this.dataset = dataset0;
		this.options = options;
		this.tdp = tdp;
		this.currentSample = new ArrayList<RootedTree>(this.nCSMC);
		for (int i = 0; i < this.nCSMC; i++)
			this.currentSample.add(initrt);

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



	public List<RootedTree> getRootedTree() {
		return currentSample;
	}
	

	public void next(Random rand) {
		iter++;
		List<RootedTree> previousSample = currentSample;
		if(sampleTrans2tranv) MHTrans2tranv(trans2tranv,  rand);

		// sample from PF
		StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();
		if ((iter % sampleTreeEveryNIter) == 0) {
			currentSample = iPmcmcStep(rand, previousSample);
			// update tdp
			if (processTree) {
				for (int k = 0; k < this.nCSMC; k++)
					tdp.process(currentSample.get(k));

			}
			if (useTopologyProcessor) {
				for (int k = 0; k < this.nCSMC; k++)
					trTopo.process(currentSample.get(k));
			}
			treeCount++;
			if (saveTreesFromPMCMC) {
				String stringOfTree = RootedTree.Util.toNewick(currentSample
						.get(0));
				String cmdStr = "echo -n 'TREE gsTree_" + treeCount + "=' >>"
						+ nameOfAllTrees;
				IO.call("bash -s", cmdStr, output);
				cmdStr = "echo '" + stringOfTree
						+ "' | sed 's/internal_[0-9]*_[0-9]*//g'"
						+ " | sed 's/leaf_//g'" + " >> " + nameOfAllTrees;
				IO.call("bash -s", cmdStr, output);
			}
		}
		// log some stats
		final int tSize = currentSample.get(0).topology().nLeaves();
		outMan.write(
				"PGS",
				"Iter",
				iter,
				"treeSize",
				tSize,
				"trans2tranv",
				trans2tranv,
				"LogLikelihood", LogLL);
	}

	private Pair<PartialCoalescentState4BackForwardKernel, Double> oneiPmcmcStep(Random rand,
			boolean isconditional,
			RootedTree onesample) {
		StoreProcessor<PartialCoalescentState4BackForwardKernel> pro = new StoreProcessor<PartialCoalescentState4BackForwardKernel>();		 
		//StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();
		CTMC currentctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), trans2tranv);	
		PartialCoalescentState init0 = PartialCoalescentState.initFastState(dataset, currentctmc, true);
		PartialCoalescentState4BackForwardKernel init = new PartialCoalescentState4BackForwardKernel(init0, init0, init0, null, 0,new int[]{-1,-1});
		//LazyParticleKernel pk2 = new BackForwardKernel(init);
		//PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, currentctmc, isGS4Clock);
		//try to modify this part, include the BF kernel.
		LazyParticleKernel kernel = new BackForwardKernel(init);
		ParticleFilter<PartialCoalescentState4BackForwardKernel> pf = new ParticleFilter<PartialCoalescentState4BackForwardKernel>();
		//ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
		pf.N = options.nParticles;
		pf.rand = rand;
		pf.resamplingStrategy = ResamplingStrategy.ALWAYS;
		//System.out.print(onesample);
		if(isconditional & (conditionedSampled!=null)){
		//List<Pair<PartialCoalescentState4BackForwardKernel, Double>> restorePCS = restoreSequence(
		//		kernel, onesample, isGS4Clock);
		List<Pair<PartialCoalescentState4BackForwardKernel, Double>> restorePCS = PartialCoalescentState4BackForwardKernel.restoreSequence(conditionedSampled);
		List<PartialCoalescentState4BackForwardKernel> path = list();
			double[] weights = new double[restorePCS.size()];
		for (int i = 0; i < restorePCS.size(); i++) {
			path.add(restorePCS.get(i).getFirst());
			weights[i] = restorePCS.get(i).getSecond();
				//System.out.print(weights[i] + "	");
		}
			//System.out.println();
		// set the conditioning and its weights
			pf.setConditional(path, weights); // here weights are in log scale.
			
			// System.out.print(restorePCS.get(restorePCS.size() -
			// 1).getFirst()); // plot
																				// the
																				// restored
																				// tree.
		}
		// do the sampling
		pf.sample(kernel, pro);
		double logMarginalLike = pf.estimateNormalizer();
		PartialCoalescentState4BackForwardKernel sampled = pro.sample(rand);
		PartialCoalescentState4BackForwardKernel conditionedSampled = sampled;
				


		return Pair.makePair(sampled, logMarginalLike);
	}


	private List<RootedTree> iPmcmcStep(Random rand,
			List<RootedTree> conditionedSamples) {
		List<RootedTree> result = new ArrayList<RootedTree>();
		List<PartialCoalescentState4BackForwardKernel> resultUnCondiPlusOne = new ArrayList<PartialCoalescentState4BackForwardKernel>(
				this.nUCSMC + 1);
		double[] logMargLikeVec = new double[this.nUCSMC + 1];
		for (int i = 0; i < this.nUCSMC; i++) {
			Pair<PartialCoalescentState4BackForwardKernel, Double> re = oneiPmcmcStep(rand,
					false, null);
			resultUnCondiPlusOne.add(re.getFirst());
			logMargLikeVec[i] = re.getSecond();
		}
		for (int i = 0; i < this.nCSMC; i++) {
			//System.out.println("Number " + i + "CSMC");
			Pair<PartialCoalescentState4BackForwardKernel, Double> re = oneiPmcmcStep(rand, true,
					conditionedSamples.get(i));
			logMargLikeVec[this.nUCSMC] = re.getSecond();
			resultUnCondiPlusOne.add(this.nUCSMC, re.getFirst());
			final double[] normalizedWeights0 = logMargLikeVec.clone();
			for (int k = 0; k < this.nUCSMC + 1; k++) {
				//System.out.print(normalizedWeights0[k] + "		");

			}
			//System.out.println();
			NumUtils.expNormalize(normalizedWeights0);
			List<Double> w = new ArrayList<Double>(this.nUCSMC + 1);
			for (int k = 0; k < this.nUCSMC + 1; k++)
 {
				w.add(k, normalizedWeights0[k]);
				//System.out.print(w.get(k) + "		");

			}
			//System.out.println();

			final int idx = Sampling.sample(rand, w);
			//System.out.println("idx is " + idx);
			result.add(resultUnCondiPlusOne.get(idx).getCurrentState().getFullCoalescentState());
		}
		return result;
	}

	private void updateLogLikelihood() {
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), trans2tranv);
		for (int i = 0; i < this.nCSMC; i++) {
		UnrootedTreeState previousncs = UnrootedTreeState.initFastState(
					currentSample.get(i).getUnrooted(), dataset, ctmc);

			previousLogLLEstimate[i] = previousncs.logLikelihood();
		}
	}

	private double sumLogLikelihood(CTMC ctmc) {
		double sumLogLikelihood = 0;
		for (int i = 0; i < this.nCSMC; i++)
			sumLogLikelihood += UnrootedTreeState.initFastState(
					currentSample.get(i).getUnrooted(), dataset, ctmc)
					.logLikelihood();
		return sumLogLikelihood;

	}

	private double sumPreviousLogLLEstimate()
	{
		double sum = 0;
		for(int i=0;i<previousLogLLEstimate.length;i++)
			sum += previousLogLLEstimate[i];
		return sum;
	}
	
	private double MHTrans2tranv(double currentTrans2tranv, Random rand) {
		Trans2tranvProposal kappaProposal=new Trans2tranvProposal(a,rand);
		Pair<Double,Double> proposed=kappaProposal.propose(currentTrans2tranv);		
		double proposedTrans2tranv=proposed.getFirst();
		RootedTree sampledTree = currentSample.get(rand.nextInt(nCSMC));
		CTMC currentctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), currentTrans2tranv);
		UnrootedTreeState ncs = UnrootedTreeState.initFastState(sampledTree.getUnrooted(), dataset, currentctmc);
		double previouslogLL = ncs.logLikelihood();
//		double previouslogLL = UnrootedTreeState.initFastState(
//				sampledTree.getUnrooted(), dataset, currentctmc)
//				.logLikelihood();
		LogLL = previouslogLL;
		
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);			
		double Loglikelihood = UnrootedTreeState.initFastState(
				sampledTree.getUnrooted(), dataset, ctmc)
				.logLikelihood();
		//UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);
		double logratio = Loglikelihood - previouslogLL+proposed.getSecond();
		//double logratio = sumLogLikelihood(ctmc) - sumPreviousLogLLEstimate() + proposed.getSecond();
		double acceptPr = Math.min(1, Math.exp(logratio));
		final boolean accept = Sampling.sampleBern(acceptPr, rand);
		if (accept) {
			trans2tranv = proposedTrans2tranv;
		//	previousLogLLEstimate = ncs.logLikelihood();			
		}
		return acceptPr;
	}


	public static double height(Map<Taxon, Double> branchLengths,
			Arbre<Taxon> arbre) {
		double sum = 0;
		do {
			arbre = arbre.getChildren().get(0);
			sum += branchLengths.get(arbre.getContents());
		} while (!arbre.isLeaf());
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
			Map.Entry entry = (Map.Entry) it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}

//	public static List<Pair<PartialCoalescentState, Double>> restoreSequence4NonClockTree(
//			PartialCoalescentState current, RootedTree rt) {
//		List<Pair<PartialCoalescentState, Double>> result = list();
//		List<Arbre<Taxon>> childrenList = rt.topology().nodes();
//		class ArbreComparator implements Comparator<Arbre<Taxon>> {
//			@Override
//			public int compare(Arbre<Taxon> arbre1, Arbre<Taxon> arbre2) {
//				return arbre1.getContents().toString()
//						.compareTo(arbre2.getContents().toString());
//			}
//		}
//		Collections.sort(childrenList, new ArbreComparator());
//		Map<Taxon, Double> branchLengths = rt.branchLengths();
//		for (int i = 0; i < childrenList.size(); i++) {
//			if (!childrenList.get(i).isLeaf()) {
//				Arbre<Taxon> currentArbre = childrenList.get(i);
//				Taxon first = currentArbre.getChildren().get(0).getContents(), second = currentArbre
//						.getChildren().get(1).getContents();
//				PartialCoalescentState coalesceResult = current.coalesce(
//						current.indexOf(first), current.indexOf(second), 0,
//						branchLengths.get(first), branchLengths.get(second),
//						currentArbre.getContents());
//				double logWeight = coalesceResult.logLikelihoodRatio(); // coalesceResult.logLikelihood()-current.logLikelihood();
//				//- Math.log(coalesceResult.nNonTrivialRoots());  
//				current = coalesceResult;
//				result.add(Pair.makePair(coalesceResult, logWeight));
//			}
//		}
//		return result;
//	}

//	public static List<Pair<PartialCoalescentState4BackForwardKernel, Double>> restoreSequence(
//			ParticleKernel<PartialCoalescentState4BackForwardKernel> kernel, RootedTree rt,
//			boolean isClock) {
//		if (!isClock)
//			throw new RuntimeException();
//			//return restoreSequence4NonClockTree(kernel.getInitial(), rt);
//		List<String> newNodeNames = list();
//		List<Pair<PartialCoalescentState4BackForwardKernel, Double>> result = list();
//		PartialCoalescentState4BackForwardKernel current = kernel.getInitial();
//		List<Arbre<Taxon>> childrenList = rt.topology().nodes();
//		Map<Taxon, Double> branchLengths = rt.branchLengths();
//		Map<Arbre<Taxon>, Double> heightMap = map();
//		for (int i = 0; i < childrenList.size(); i++) {
//			if (!childrenList.get(i).isLeaf())
//				heightMap.put(childrenList.get(i),
//						height(branchLengths, childrenList.get(i)));
//		}
//		// sort the height in an increasing order
//		Map<Arbre<Taxon>, Double> sortHeightMap = sortByValue(heightMap);
//		Set<Arbre<Taxon>> arbreSet = sortHeightMap.keySet();
//		Iterator<Arbre<Taxon>> arbreIterator = arbreSet.iterator();
//		while (arbreIterator.hasNext()) {
//			Arbre<Taxon> currentArbre = arbreIterator.next();
//			newNodeNames.add(currentArbre.toString());
//		}
//		arbreIterator = arbreSet.iterator();
//		double previousHeight = 0;
//		while (arbreIterator.hasNext()) {
//			Arbre<Taxon> currentArbre = arbreIterator.next();
//			Taxon first = currentArbre.getChildren().get(0).getContents(), second = currentArbre
//					.getChildren().get(1).getContents();
//			double currentHeight = sortHeightMap.get(currentArbre);
//			double currentDelta = currentHeight - previousHeight;
//			previousHeight = currentHeight;
//			PartialCoalescentState coalesceResult = current.getCurrentState().coalesce(
//					current.getCurrentState().indexOf(first), current.getCurrentState().indexOf(second),
//					currentDelta, 0, 0, currentArbre.getContents());
//			double logWeight = coalesceResult.logLikelihood()
//					- current.getCurrentState().logLikelihood();
//
//			// coalesceResult.logLikelihoodRatio(); // equal to
//																	// coalesceResult.logLikelihood()
//																	// -
//																	// current.logLikelihood();
//			current = coalesceResult;
//			result.add(Pair.makePair(coalesceResult, logWeight));
//		}
//		return result;
//	}
//
//
//
}
