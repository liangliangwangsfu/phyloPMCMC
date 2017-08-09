package phyloPMCMC;

import java.util.*;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.io.Dataset.DatasetUtils;
import pty.io.TreeEvaluator.TreeMetric;
import pty.smc.LazyPCS;
import pty.smc.LazyParticleFilter;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.models.CTMC;


import ev.poi.processors.TreeDistancesProcessor;
import gep.util.OutputManager;
import ma.SequenceType;
import nuts.math.Sampling;
import nuts.util.CollUtils;


public class PMMHSMC2
{
	private final Dataset dataset;
	ParticleFilterOptions options=null;
	private final TreeDistancesProcessor tdp;
	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
	private RootedTree currentSample = null;
	private double trans2tranv;
	public double a=1.2;
	public static OutputManager outMan = new OutputManager();	
	public boolean stop=false;	
	private int iter=0;
	public double priorRate=2.0; 
	private int tryCountBeforAccept=0;


	public PMMHSMC2(Dataset dataset0, ParticleFilterOptions options,  TreeDistancesProcessor tdp, double trans2tranv0)
	{
		this.options=options;
		this.trans2tranv=trans2tranv0;
		this.dataset=dataset0;
		this.tdp = tdp;
	}

	
	public void next(Random rand)
	{
		iter++;
		//		System.out.println(rand.nextDouble());
		RootedTree previousSample = currentSample;		
		double scale=Sampling.nextDouble(rand, 1.0/a, a);
		//		System.out.println(scale);
		double proposedTrans2tranv = scale*trans2tranv; 
		//		double proposedTrans2tranv = 2.0;
		// sample from PF
		StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();		 
//		TreeDistancesProcessor innerTdp = new TreeDistancesProcessor();
		
		double acceptPr = 0.0;
		try 
		{
			tryCountBeforAccept++;
//			if(tryCountBeforAccept>10) stop=true;

			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);  
			PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, ctmc, true);  // is clock	          								
			LazyParticleKernel kernel = new PriorPriorKernel(init);
			LazyParticleFilter<PartialCoalescentState> pf = new LazyParticleFilter<PartialCoalescentState>(kernel, options);		    		
			double zHat=pf.sample(pro);
			
			// compute accept/reject
			double marginalLoglike=zHat+init.logLikelihood();
			final double logRatio =marginalLoglike  - previousLogLLEstimate+ trans2tranv*(1-scale)/priorRate + Math.log(scale);
			acceptPr = Math.min(1, Math.exp(logRatio));
			if (currentSample != null &&  Double.isInfinite(acceptPr))
				throw new RuntimeException();
		
			final boolean accept = Sampling.sampleBern(acceptPr, rand);
			if (accept)
			{
				// sample from sample
				PartialCoalescentState sampled = pro.sample(rand);				
				currentSample = sampled.getFullCoalescentState();
				// set ll!
				previousLogLLEstimate = marginalLoglike;
				trans2tranv=proposedTrans2tranv;
				
			}
		}
		catch (Exception e)
		{
			// total fail!
		}
		// update tdp
		tdp.process(currentSample);
		// log some stats
		final int tSize = currentSample.topology().nLeaves();
		outMan.write("PMMH",
				"Iter", iter,
				"treeSize", tSize,
				"acceptPr", acceptPr, 
				//        "maskSparsity", currentSparsity,
				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
				"trans2tranv", trans2tranv, 
				"LogLikelihood", previousLogLLEstimate);
	}
}
